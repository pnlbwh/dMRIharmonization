#!/usr/bin/env python

from plumbum import cli
from distutils.spawn import find_executable
import os, shutil
from dti import dti
from rish import rish
from buildTemplate import difference_calc, antsMult, warp_bands, \
    dti_stat, rish_stat, template_masking, createAntsCaselist, \
    read_caselist
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
from debug_fa import sub2tmp2mni, analyzeStat
import numpy as np

def check_csv(file, force):

    with open(file) as f:
        content= f.read()

        for line, row in enumerate(content.split()):
            dwi_mask= [element for element in row.split(',') if element] # handling w/space
            if len(dwi_mask) != 2:
                raise FileNotFoundError(f'Columns don\'t have same number of entries: check line {line} in {file}')

            dirCheckFlag= 1
            for img in dwi_mask:
                if not os.path.exists(img):
                    raise FileNotFoundError(f'{img} does not exist: check line {line} in {file}')

                elif dirCheckFlag:
                    # create DTI and harmonization directory
                    dtiPath= os.path.join(os.path.dirname(img),'dti')
                    check_dir(dtiPath, force)

                    harmPath= os.path.join(os.path.dirname(img),'harm')
                    check_dir(harmPath, force)

                    # rishPath= os.path.join(os.path.dirname(img), 'harm', 'rish')
                    # check_dir(rishPath, force)

                    dirCheckFlag= 0


def check_dir(path, force):
    if os.path.exists(path) and force:
        warnings.warn(f'{path} exists and will be overwritten')
        shutil.rmtree(path)
        os.makedirs(path)
    elif not os.path.exists(path):
        os.makedirs(path)
    else:
        warnings.warn(f'{path} exists, --force not specified, continuing with existing directory')



def dti_harm(imgPath, maskPath, N_shm):

    directory = os.path.dirname(imgPath)
    inPrefix = imgPath.split('.')[0]
    prefix = os.path.split(inPrefix)[-1]

    outPrefix = os.path.join(directory, 'dti', prefix)

    # if the dti output exists with the same prefix, don't dtifit again
    if not os.path.exists(outPrefix+'_FA.nii.gz'):
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

    harm_csv = cli.SwitchAttr(
        ['--harmonized'],
        cli.ExistingFile,
        help='harmonized csv/txt file with first column for dwi and 2nd column for mask: dwi1,mask1\ndwi2,mask2\n...',
        mandatory=False)

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

    debug = cli.Flag(
        '--debug',
        help= 'turn on this flag to debug harmonized data (valid only with --process)',
        default= False)

    reference = cli.SwitchAttr(
        '--refName',
        help= 'reference site name',
        mandatory= True)

    target = cli.SwitchAttr(
        '--targetName',
        help= 'target site name',
        mandatory= True)

    diffusionMeasures = ['MD', 'FA', 'GFA']

    def preprocessing(self, imgPath, maskPath):

        # load signal attributes for pre-processing ----------------------------------------------------------------
        lowRes = load(imgPath)
        lowResImg = lowRes.get_data()
        lowResImgHdr = lowRes.header

        lowRes = load(maskPath)
        lowResMask = lowRes.get_data()
        lowResMaskHdr = lowRes.header

        lowResImg = applymask(lowResImg, lowResMask)

        inPrefix = imgPath.split('.')[0]

        bvals, _ = read_bvals_bvecs(inPrefix + '.bval', None)

        # pre-processing -------------------------------------------------------------------------------------------
        suffix= None
        # modifies data only
        if self.denoise:
            print('Denoising ', imgPath)
            lowResImg, _ = denoising(lowResImg, lowResMask)
            suffix = '_denoised'
            # save_nifti(imgPath.split('.')[0]+'_denoised.nii.gz', lowResImg, lowResImgHdr.affine)

        # modifies data, and bvals
        newBval= float(self.bvalMap)
        if self.bvalMap and (max(bvals)>1.01*newBval or max(bvals)<0.99*newBval):
            print('B value mapping ', imgPath)
            lowResImg, bvals = bvalMap(lowResImg, bvals, newBval)
            suffix = '_bmapped'
            # save_nifti(imgPath.split('.')[0]+'_bmapped.nii.gz', lowResImg, lowResImgHdr.affine)

        # modifies data, mask, and headers
        if self.resample:
            print('Resampling ', imgPath)
            sp_high = np.array([float(i) for i in self.resample.split('x')])
            imgPath, maskPath = \
                resampling(imgPath, maskPath, lowResImg, lowResImgHdr, lowResMask, lowResMaskHdr, sp_high, bvals)
            suffix = '_resampled'

        # save pre-processed data; resampled data is saved inside resampling() -------------------------------------
        if (self.denoise or self.bvalMap) and not self.resample:
            imgPath = inPrefix + suffix + '.nii.gz'
            save_nifti(imgPath, lowResImg, lowResImgHdr.get_qform())

        if suffix:
            shutil.copyfile(inPrefix + '.bvec', inPrefix + suffix + '.bvec')
        if self.bvalMap:
            write_bvals(inPrefix + suffix + '.bval', bvals)
        elif self.denoise or self.resample:
            shutil.copyfile(inPrefix + '.bval', inPrefix + suffix + '.bval')


        return (imgPath, maskPath)


    def common_processing(self, caselist):

        imgs, masks = read_caselist(caselist)
        f= open(caselist+'.modified', 'w')

        # TODO: parellelize
        for i in range(len(imgs)):

            imgs[i], masks[i]= self.preprocessing(imgs[i], masks[i])
            dti_harm(imgs[i], masks[i], self.N_shm)
            f.write(f'{imgs[i]},{masks[i]}\n')

        f.close()

        return (imgs, masks)



    def createTemplate(self):

        # check directory existence
        check_dir(self.templatePath, self.force)

        # go through each file listed in csv, check their existence, create dti and harm directories
        check_csv(self.ref_csv, self.force)
        check_csv(self.target_csv, self.force)

        # createTemplate steps -----------------------------------------------------------------------------------------

        # read image lists
        refImgs, refMasks= self.common_processing(self.ref_csv)
        if not self.ref_csv.endswith('.modified'):
            self.ref_csv += '.modified'
        # debug: use the following line to omit processing again
        # refImgs, refMasks = read_caselist(self.ref_csv)

        targetImgs, targetMasks= self.common_processing(self.target_csv)
        if not self.target_csv.endswith('.modified'):
            self.target_csv += '.modified'
        # debug: use the following line to omit processing again
        # targetImgs, targetMasks = read_caselist(self.target_csv)

        imgs= refImgs+targetImgs
        masks= refMasks+targetMasks

        # create caselist for antsMult
        antsMultCaselist= os.path.join(self.templatePath, 'antsMultCaselist.txt')
        createAntsCaselist(imgs, antsMultCaselist)
        # run ANTS multivariate template construction
        antsMult(antsMultCaselist, self.templatePath)

        # load templateAffine
        templateAffine= load_nifti(os.path.join(self.templatePath, 'template0.nii.gz'))[1]

        # warp mask, dti, and rish bands
        # TODO: parellelize
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
        check_csv(self.target_csv, self.force)

        # cleanOutliers steps ------------------------------------------------------------------------------------------

        # read target image list
        moving= os.path.join(self.templatePath, f'Mean_{self.target}_FA.nii.gz')
        imgs, masks= read_caselist(self.target_csv)

        preFlag= 1 # omit preprocessing of target data again
        if self.target_csv.endswith('.modified'):
            preFlag= 0
        else:
            # this file will be used later for debugging
            self.target_csv += '.modified'
            fm = open(self.target_csv, 'w')


        # TODO: parellelize
        self.harm_csv= self.target_csv+'.harmonized'
        fh= open(self.harm_csv, 'w')
        for imgPath, maskPath in zip(imgs, masks):

            if preFlag:
                imgPath, maskPath = self.preprocessing(imgPath, maskPath)
                fm.write(imgPath + ',' + maskPath + '\n')

            b0, shm_coeff, qb_model = dti_harm(imgPath, maskPath, self.N_shm)

            directory= os.path.dirname(imgPath)
            inPrefix= imgPath.split('.')[0]
            prefix= os.path.split(inPrefix)[-1]

            print(f'Registering template FA to {imgPath} space ...')
            outPrefix= os.path.join(directory, 'harm', 'ToSubjectSpace_'+ prefix)
            fixed= os.path.join(directory, 'dti', f'{prefix}_FA.nii.gz')
            antsReg(fixed, maskPath, moving, outPrefix)
            antsApply(self.templatePath, os.path.join(directory, 'harm'), prefix, self.N_shm)

            print(f'Reconstructing signal from {imgPath} rish features ...')
            harmImg, harmMask= ring_masking(directory, prefix, maskPath, self.N_shm, shm_coeff, b0, qb_model)
            fh.write(harmImg+','+harmMask+'\n')

            shutil.copyfile(inPrefix + '.bvec', harmImg.split('.')[0] + '.bvec')
            shutil.copyfile(inPrefix + '.bval', harmImg.split('.')[0] + '.bval')

            if self.debug:
                dti_harm(harmImg, harmMask, self.N_shm)

        if preFlag:
            fm.close()
        fh.close()
        print('\n\nHarmonization completed\n\n')


    def post_debug(self):

        print('\n\n Reference site')
        sub2tmp2mni(self.templatePath, self.reference, self.ref_csv, self.diffusionMeasures)
        ref_mean = analyzeStat(self.ref_csv)

        print('\n\n Target site before harmonization')
        sub2tmp2mni(self.templatePath, self.target, self.target_csv, self.diffusionMeasures)
        target_mean_before = analyzeStat(self.target_csv)

        print('\n\n Target site after harmonization')
        sub2tmp2mni(self.templatePath, self.target, self.harm_csv, self.diffusionMeasures)
        target_mean_after = analyzeStat(self.harm_csv)

        print('\n\nPrinting statistics :')
        print(f'{self.reference} mean FA: ', np.mean(ref_mean))
        print(f'{self.target} mean FA before harmonization: ', np.mean(target_mean_before))
        print(f'{self.target} mean FA after harmonization: ', np.mean(target_mean_after))


    def sanityCheck(self):

        if not (self.create or self.process or self.debug):
            raise AttributeError('No option selected, ' 
                                'specify one (or many of) creation, harmonization, and debug flags')

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

        if self.debug:
            self.post_debug()


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


Harmonization tests:


python -m pdb \
/home/tb571/Downloads/Harmonization-Python/lib/harmonization.py \
--N_shm 6 \
--resample 1.5x1.5x1.5 \
--bvalMap 1000 \
--target /home/tb571/Downloads/Harmonization-Python/BSNIP_Baltimore/target_caselist.txt.modified \
--reference /home/tb571/Downloads/Harmonization-Python/BSNIP_Baltimore/ref_caselist.txt.modified \
--refName CIDAR \
--targetName BSNIP \
--template /home/tb571/Downloads/Harmonization-Python/BSNIP_Baltimore/template/ \
--create \
--process \
--debug


python -m pdb \
/home/tb571/Downloads/Harmonization-Python/lib/harmonization.py \
--N_shm 6 \
--target /home/tb571/Downloads/Harmonization-Python/connectom_prisma_demoData/target_caselist.txt \
--reference /home/tb571/Downloads/Harmonization-Python/connectom_prisma_demoData/ref_caselist.txt \
--refName CIDAR \
--targetName BSNIP \
--template /home/tb571/Downloads/Harmonization-Python/connectom_prisma_demoData/template/ \
--create \
--process \
--debug

'''
