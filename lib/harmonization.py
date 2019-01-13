#!/usr/bin/env python

from plumbum import cli
from distutils.spawn import find_executable
import os, shutil
from dti import dti
from rish import rish
from buildTemplate import difference_calc, antsMult, warp_bands, \
    dti_stat, rish_stat, template_masking, createAntsCaselist
from denoising import denoising
from bvalMap import bvalMap
from resampling import resampling
import warnings
with warnings.catch_warnings():
    warnings.filterwarnings("ignore", category=FutureWarning)

    from dipy.io.image import load_nifti, save_nifti
    from dipy.io import read_bvals_bvecs
    from dipy.segment.mask import applymask
    from nibabel import load

from cleanOutliers import antsReg, antsApply, ring_masking
import numpy as np

def check_csv(file, force):

    with open(file) as f:
        content= f.read()

        for line, row in enumerate(content.split()):
            dwi_mask= [element for element in row.split(',') if element] # handling w/space
            if len(dwi_mask) != 2:
                raise FileNotFoundError(f'Columns don\'t have same number of entries: check line {line} in {file}')

            for img in dwi_mask:
                if not os.path.exists(img):
                    raise FileNotFoundError(f'{img} does not exist: check line {line} in {file}')

                else:
                    # create DTI and harmonization directory
                    dtiPath= os.path.join(os.path.dirname(img),'dti')
                    check_dir(dtiPath, force)

                    harmPath= os.path.join(os.path.dirname(img),'harm')
                    check_dir(harmPath, force)

                    # rishPath= os.path.join(os.path.dirname(img), 'harm', 'rish')
                    # check_dir(rishPath, force)


def check_dir(path, force):
    if os.path.exists(path) and force:
        warnings.warn(f'{path} exists and will be overwritten')
        shutil.rmtree(path)
        os.makedirs(path)
    elif not os.path.exists(path):
        os.makedirs(path)
    else:
        raise IsADirectoryError(f'{path} exists, use --force to overwrite or delete existing directory')


def read_caselist(file):

    with open(file) as f:

        imgs = []
        masks = []
        content= f.read()
        for line, row in enumerate(content.split()):
            temp= [element for element in row.split(',') if element] # handling w/space
            imgs.append(temp[0])
            masks.append(temp[1])

        return (imgs, masks)


# def template_csv(file1, file2, templatePath):
#
#     f1 = open(file1)
#     f2 = open(file2)
#     file3= os.path.join(templatePath,'ref_target_caselist.txt')
#     f3 = open(file3, 'w')
#
#     content = f1.read()+f2.read()
#
#     f3.write(content)
#
#     f1.close()
#     f2.close()
#     f3.close()
#
#     return file3


def dti_harm(imgPath, maskPath, N_shm):

    directory = os.path.dirname(imgPath)
    inPrefix = imgPath.split('.')[0]
    prefix = os.path.split(inPrefix)[-1]

    outPrefix = os.path.join(directory, 'dti', prefix)
    dti(imgPath, maskPath, inPrefix, outPrefix)

    outPrefix = os.path.join(directory, 'harm', prefix)
    b0, shm_coeff, qb_model= rish(imgPath, maskPath, inPrefix, outPrefix, N_shm)

    return (b0, shm_coeff, qb_model)


def write_bvals(bval_file, bvals):
    with open(bval_file, 'w') as f:
        f.write(('\n').join(str(b) for b in bvals))


class pipeline(cli.Application):

    """Template creation and harmonization"""

    ref_csv = cli.SwitchAttr(
        ['--reference'],
        cli.ExistingFile,
        help='reference csv/txt file with first column for dwi and 2nd column for mask: dwi1,mask1\ndwi2,mask2\n...',
        mandatory=False)

    target_csv = cli.SwitchAttr(
        ['--target'],
        cli.ExistingFile,
        help='target csv/txt file with first column for dwi and 2nd column for mask: dwi1,mask1\ndwi2,mask2\n...',
        mandatory=True)

    templatePath = cli.SwitchAttr(
        ['--template'],
        help='template directory',
        mandatory=True)

    # dtiPath = cli.SwitchAttr(
    #     ['--dti'],
    #     help='dti directory',
    #     mandatory=True)
    #
    # harmPath = cli.SwitchAttr(
    #     ['--harm'],
    #     help='harmonization directory',
    #     mandatory=True)

    N_shm = cli.SwitchAttr(
        ['--N_shm'],
        help='spherical harmonic order',
        default= 8)

    force = cli.Flag(
        ['--force'],
        help='turn on this flag to overwrite existing data',
        default= False)

    travelHeads = cli.Flag(
        ['--travelHeads'],
        help='travelling heads',
        default= False)

    resample = cli.SwitchAttr(
        '--resample',
        help='voxel size MxNxO to resample into',
        default= False)

    bvalMap = cli.SwitchAttr(
        '--bvalMap',
        help='specify a bmax to scale bvalues into',
        default= False)

    denoise = cli.Flag(
        '--denoise',
        help='turn on this flag to denoise voxel data',
        default= False)

    create = cli.Flag(
        '--create',
        help= 'turn on this flag to create template',
        default= False)

    process = cli.Flag(
        '--process',
        help= 'turn on this flag to harmonize',
        default= False)

    reference= 'reference'
    target= 'target'
    diffusionMeasures = ['MD', 'FA']

    def preprocessing(self, imgPath, maskPath):

        # load signal attributes for pre-processing ----------------------------------------------------------------
        lowRes = load(imgPath)
        lowResImg = lowRes.get_data()
        lowResImgHdr = lowRes.header

        lowRes = load(maskPath)
        lowResMask = lowRes.get_data()
        lowResMaskHdr = lowRes.header

        lowResImg = applymask(lowResImg, lowResMask)

        # directory = os.path.dirname(imgPath)
        inPrefix = imgPath.split('.')[0]
        # prefix = os.path.split(inPrefix)[-1]

        bvals, _ = read_bvals_bvecs(inPrefix + '.bval', None)

        # pre-processing -------------------------------------------------------------------------------------------

        # modifies data only
        suffix = '_denoised'
        imgPath = inPrefix + suffix + '.nii.gz'
        fileExist= os.path.exists(imgPath)
        if self.denoise and (not fileExist or (fileExist and self.force)):
            print('Denoising ', imgPath)
            lowResImg, _ = denoising(lowResImg, lowResMask)
            # suffix = '_denoised'
            # save_nifti(imgPath.split('.')[0]+'_denoised.nii.gz', lowResImg, lowResImgHdr.affine)
        elif fileExist and not self.force:
            raise FileExistsError(f'Denoised {imgPath} exists, use --force to overwrite or delete existing directory')

        # modifies data, and bvals
        suffix = '_bmapped'
        imgPath = inPrefix + suffix + '.nii.gz'
        fileExist= os.path.exists(imgPath)
        if self.bvalMap and (not fileExist or (fileExist and self.force)):
            print('B value mapping ', imgPath)
            lowResImg, bvals = bvalMap(lowResImg, bvals, float(self.bvalMap))
            # suffix = '_bmapped'
            # write_bvals(inPrefix + suffix + '.bval', bvals)
            # save_nifti(imgPath.split('.')[0]+'_bmapped.nii.gz', lowResImg, lowResImgHdr.affine)
        elif fileExist and not self.force:
            raise FileExistsError(f'B value mapped {imgPath} exists, use --force to overwrite or delete existing directory')

        # modifies data, mask, and headers
        suffix = '_resampled'
        imgPath = inPrefix + suffix + '.nii.gz'
        fileExist= os.path.exists(imgPath)
        if self.resample and (not fileExist or (fileExist and self.force)):
            print('Resampling ', imgPath)
            sp_high = np.array([float(i) for i in self.resample.split('x')])
            imgPath, maskPath = \
                resampling(imgPath, maskPath, lowResImg, lowResImgHdr, lowResMask, lowResMaskHdr, sp_high, bvals)
            # suffix = '_resampled'
        elif fileExist and not self.force:
            raise FileExistsError(f'Resampled {imgPath} exists, use --force to overwrite or delete existing directory')

        # save pre-processed data; resampled data is saved inside resampling() -------------------------------------
        if (self.denoise or self.bvalMap) and not self.resample:
            imgPath = inPrefix + suffix + '.nii.gz'
            save_nifti(imgPath, lowResImg, lowResImgHdr.affine)

        shutil.copyfile(inPrefix + '.bvec', inPrefix + suffix + '.bvec')
        if self.bvalMap:
            write_bvals(inPrefix + suffix + '.bval', bvals)
        elif self.denoise or self.resample:
            shutil.copyfile(inPrefix + '.bval', inPrefix + suffix + '.bval')


        return (imgPath, maskPath)


    def common_processing(self, caselist):

        imgs, masks = read_caselist(caselist)
        f= open(caselist+'.modified', 'w')

        for i in range(len(imgs)):

            imgs[i], masks[i]= self.preprocessing(imgs[i], masks[i])

            dti_harm(imgs[i], masks[i], self.N_shm)

            f.write(f'{imgs[i]},{masks[i]}\n')

        f.close()

        return (imgs, masks)



    def createTemplate(self):

        # go through each file listed in csv, check their existence, create dti and harm directories
        # if self.ref_csv:
        #     check_csv(self.ref_csv, self.force)

        # check directory existence
        # check_dir(self.templatePath, self.force)

        # dtifit and rish feature
        # imgs, masks= self.common_processing(allCaselist)
        # debug: use the following line to omit processing again
        # imgs, masks= read_caselist(self.ref_csv)

        # create extended caselist
        # allCaselist= template_csv(self.ref_csv, self.target_csv, self.templatePath)

        # createTemplate steps -----------------------------------------------------------------------------------------

        # read image lists
        # refImgs, refMasks= self.common_processing(self.ref_csv)
        self.ref_csv+='.modified'
        refImgs, refMasks = read_caselist(self.ref_csv)

        # targetImgs, targetMasks= self.common_processing(self.target_csv)
        self.target_csv+='.modified'
        targetImgs, targetMasks = read_caselist(self.target_csv)
        imgs= refImgs+targetImgs
        masks= refMasks+targetMasks

        # create caselist for antsMult
        # antsMultCaselist= os.path.join(self.templatePath, 'antsMultCaselist.txt')
        # createAntsCaselist(imgs, antsMultCaselist)
        # # run ANTS multivariate template construction
        # antsMult(antsMultCaselist, self.templatePath)

        # load templateAffine
        templateAffine= load_nifti(os.path.join(self.templatePath, 'template0.nii.gz'))[1]

        # warp mask, dti, and rish bands
        for imgPath, maskPath in zip(imgs, masks):
            warp_bands(imgPath, maskPath, self.templatePath, self.N_shm, self.diffusionMeasures)


        print('dti statistics: mean, std(FA, MD) calculation of reference site')
        refMaskPath= dti_stat(self.reference, refImgs, refMasks, self.templatePath, templateAffine, self.diffusionMeasures)
        print('dti statistics: mean, std(FA, MD) calculation of target site')
        targetMaskPath= dti_stat(self.target, targetImgs, targetMasks, self.templatePath, templateAffine, self.diffusionMeasures)

        print('masking dti statistics of reference site')
        _= template_masking(refMaskPath, targetMaskPath, self.templatePath, self.reference, self.diffusionMeasures)
        print('masking dti statistics of target site')
        templateMask= template_masking(refMaskPath, targetMaskPath, self.templatePath, self.target, self.diffusionMeasures)

        print('rish_statistics mean, std(L{i}) calculation or reference site')
        rish_stat(self.reference, imgs, self.templatePath, templateAffine, self.N_shm)
        print('rish_statistics mean, std(L{i}) calculation or target site')
        rish_stat(self.target, imgs, self.templatePath, templateAffine, self.N_shm)

        print('difference calculation of diffusionMeasures')
        difference_calc(self.reference, self.target, refImgs, targetImgs, self.templatePath, templateAffine,
                        'dti', templateMask, self.diffusionMeasures, self.travelHeads)

        print('difference calculation of rishFeatures')
        difference_calc(self.reference, self.target, refImgs, targetImgs, self.templatePath, templateAffine,
                        'harm', templateMask, [f'L{i}' for i in range(0, self.N_shm+1, 2)], self.travelHeads)


        print('\n\nTemplate creation completed \n\n')


    def harmonizeData(self):

        # check the templatePath
        if not os.path.exists(self.templatePath):
            raise NotADirectoryError(f'{self.templatePath} does not exist')
        else:
            if not os.listdir(self.templatePath):
                raise ValueError(f'{self.templatePath} is empty')

        # go through each file listed in csv, check their existence, create dti and harm directories
        # check_csv(self.target_csv, self.force)

        # self.common_processing(self.target_csv)

        # cleanOutliers steps ------------------------------------------------------------------------------------------

        # read target image list
        self.target_csv += '.modified'
        moving= os.path.join(self.templatePath, f'Mean_{self.target}_FA.nii.gz')
        imgs, masks= read_caselist(self.target_csv)
        for imgPath, maskPath in zip(imgs, masks):

            directory= os.path.dirname(imgPath)
            inPrefix= imgPath.split('.')[0]
            prefix= os.path.split(inPrefix)[-1]

            b0, shm_coeff, qb_model = dti_harm(imgPath, maskPath, self.N_shm)

            # if self.force or not os.path.exists(inPrefix+'.npz'):
            #     imgPath, maskPath= self.preprocessing(imgPath, maskPath)
            #     b0, shm_coeff, qb_model= dti_harm(imgPath, maskPath, self.N_shm)
            # else:
            #     content= np.load(inPrefix+'.npz')
            #     b0= content['b0']
            #     shm_coeff= content['shm_coeff']
            #     qb_model= content['qb_model']

            print(f'Registering {imgPath} FA to subject space ...')
            outPrefix= os.path.join(directory, 'harm', 'ToSubjectSpace_'+ prefix)
            fixed= os.path.join(directory, 'dti', f'{prefix}_FA.nii.gz')
            antsReg(fixed, maskPath, moving, outPrefix)
            antsApply(self.templatePath, os.path.join(directory, 'harm'), prefix, self.N_shm)

            print(f'Harmonizing {imgPath} rish features ...')
            mappedFile= ring_masking(directory, prefix, maskPath, self.N_shm, shm_coeff, b0, qb_model)
            outPrefix= os.path.join(directory, 'harm', f'harmonized_{prefix}')
            rish(mappedFile, maskPath, inPrefix, outPrefix, self.N_shm, qb_model)

        print('\n\nHarmonization completed\n\n')


    def sanityCheck(self):

        if not (self.create or self.process):
            raise AttributeError('No option selected, ' 
                                'specify either creation and/or harmonization flag')

        # check ants commands
        external_commands= [
            'antsMultivariateTemplateConstruction2.sh',
            'antsApplyTransforms',
            'antsRegistrationSyNQuick.sh']

        for cmd in external_commands:
            exe= find_executable(cmd)
            if not exe:
                raise EnvironmentError(f'{cmd} not found')


        # go through each file listed in csv, check their existence, create dti and harm directories
        # if self.ref_csv:
        #     check_csv(self.ref_csv, self.force)
        # check_csv(self.target_csv, self.force)



    def main(self):

        self.N_shm= int(self.N_shm)
        self.sanityCheck()

        if self.create:
            self.createTemplate()

        if self.process:
            self.harmonizeData()


if __name__ == '__main__':
    pipeline.run()

'''

python -m pdb \
/home/tb571/Downloads/Harmonization-Python/lib/harmonization.py \
--N_shm 6 \
--denoise \
--resample 1.5x1.5x1.5 \
--bvalMap 1000 \
--force \
--create \
--target /home/tb571/Downloads/Harmonization-Python/test_data/target_caselist.txt \
--reference /home/tb571/Downloads/Harmonization-Python/test_data/ref_caselist.txt \
--template /home/tb571/Downloads/Harmonization-Python/test_data/template/ \
--travelHeads


python -m pdb \
/home/tb571/Downloads/Harmonization-Python/lib/harmonization.py \
--N_shm 6 \
--denoise \
--resample 1.5x1.5x1.5 \
--bvalMap 1000 \
--process \
--target /home/tb571/Downloads/Harmonization-Python/test_data/target_caselist.txt \
--template /home/tb571/Downloads/Harmonization-Python/test_data/template/ \
--travelHeads

'''
