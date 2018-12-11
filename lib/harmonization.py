#!/usr/bin/env python

from plumbum import cli
from distutils.spawn import find_executable
import os, warnings, shutil
from dti import dti
from rish import rish
from buildTemplate import difference_calc, antsMult, warp_bands, dti_stat, rish_stat, template_masking
import warnings
with warnings.catch_warnings():
    warnings.filterwarnings("ignore", category=FutureWarning)

    from dipy.io.image import load_nifti, save_nifti
    from dipy.reconst.shm import normalize_data

from cleanOutliers import antsReg, antsApply, ring_masking

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

                    rishPath= os.path.join(os.path.dirname(img),'harm', 'rish')
                    check_dir(rishPath, force)



def read_caselist(file):

    with open(file) as f:
        content= f.read()

        imgs= [], masks= []
        for line, row in enumerate(content.split()):
            temp= [element for element in row.split(',') if element] # handling w/space
            imgs.append(temp[0])
            masks.append(temp[1])

    return imgs, masks



def template_csv(file1, file2, templatePath):

    f1 = open(file1)
    f2 = open(file2)
    f3 = open(os.path.join(templatePath,'antsMultCaselist.txt'), 'w')

    content = f1.read()+f2.read()

    f3.write(content)

    f1.close()
    f2.close()
    f3.close()


def check_dir(path, force):

    if os.path.exists(path) and force:
        warnings.warn(f'{path} exists and will be overwritten')
        shutil.rmtree(path)
        os.makedirs(path)
    elif not os.path.exists(path):
        os.makedirs(path)
    else:
        raise IsADirectoryError(f'{path} exists, use --force to overwrite or delete existing directory')

def dti_harm(imgPath, maskPath, N_shm):

    directory = os.path.dirname(imgPath)
    inPrefix = imgPath.split('.')[0]
    prefix = os.path.split(inPrefix)[-1]

    outPrefix = os.path.join(directory, 'dti', prefix)
    dti(imgPath, maskPath, inPrefix, outPrefix)

    outPrefix = os.path.join(directory, 'harm', prefix)
    b0, shm_coeff, fit_matrix= rish(imgPath, maskPath, inPrefix, outPrefix, N_shm)

    return (b0, shm_coeff, fit_matrix)


class pipeline(cli.Application):

    """Template creation and harmonization"""

    ref_csv = cli.SwitchAttr(
        ['--reference'],
        cli.ExistingFile,
        help='csv/txt file with first column for dwi and 2nd column for mask: dwi1,mask1\ndwi2,mask2\n...',
        mandatory=False)

    target_csv = cli.SwitchAttr(
        ['--target'],
        cli.ExistingFile,
        help='csv/txt file with first column for dwi and 2nd column for mask: dwi1,mask1\ndwi2,mask2\n...',
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

    resample = cli.Flag(
        '--resample',
        help='turn on this flag to resample voxel data',
        default= False)

    bvalMap = cli.Flag(
        '--bvalMap',
        help='turn on this flag to remap the bvalues',
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

    def common_processing(self, caselist):

        imgs, masks = read_caselist(caselist)

        # preprocessing
        # dti
        # rish
        # dti

        if self.denoise:
            self.denoiseImg()

        if self.bvalMap:
            self.mapBvals()

        if self.resample():
            self.resampleImg()


        for imgPath, maskPath in (imgs, masks):
            dti_harm(imgPath, maskPath, self.N_shm)

        return (imgs, masks)

    def createTemplate(self):

        # check directory existence
        check_dir(self.templatePath, self.force)

        # create extended caselist
        template_csv(self.ref_csv, self.target_csv, self.templatePath)

        imgs, masks= self.common_processing(os.path.join(self.templatePath, 'antsMultCaselist.txt'))

        # createTemplate steps

        # read image lists
        refImgs, refMasks= read_caselist(self.ref_csv)
        targetImgs, targetMasks= read_caselist(self.target_csv)

        # run ANTS multivariate template construction
        antsMult(os.path.join(self.templatePath, 'antsMultCaselist.txt'), self.templatePath)


        # load templateAffine
        templateAffine= load_nifti(os.path.join(self.templatePath, 'template0.nii.gz'))[1]

        # warp mask, dti, and rish bands
        for imgPath, maskPath in (imgs, masks):
            warp_bands(imgPath, maskPath, self.templatePath, self.N_shm, self.diffusionMeasures)


        # dti statistics mean, std(FA, MD) calculation
        refMaskPath= dti_stat(self.reference, refImgs, refMasks, self.templatePath, templateAffine, self.diffusionMeasures)
        targetMaskPath= dti_stat(self.target, targetImgs, targetMasks, self.templatePath, templateAffine, self.diffusionMeasures)

        # masking dti statistics
        _= template_masking(refMaskPath, targetMaskPath, self.templatePath, self.reference, self.diffusionMeasures)
        templateMask= template_masking(refMaskPath, targetMaskPath, self.templatePath, self.target, self.diffusionMeasures)

        # rish_statistics mean, std(L{i}) calculation
        rish_stat(self.reference, imgs, self.templatePath, templateAffine, self.N_shm)
        rish_stat(self.target, imgs, self.templatePath, templateAffine, self.N_shm)

        # difference calculation
        difference_calc(self.reference, self.target, refImgs, targetImgs, self.templatePath, templateAffine,
                        'dti', templateMask, self.diffusionMeasures, self.travelHeads)

        difference_calc(self.reference, self.target, refImgs, targetImgs, self.templatePath, templateAffine,
                        'harm', templateMask, [f'L{i}' for i in range(0, self.N_shm, 2)], self.travelHeads)




    def harmonizeData(self):

        # check the templatePath
        if os.path.exists(self.templatePath):
            raise NotADirectoryError(f'{self.templatePath} does not exist')
        else:
            if not os.listdir(self.templatePath):
                raise ValueError(f'{self.templatePath} is empty')


        # self.common_processing(self.target_csv)

        # cleanOutliers steps

        # read target image list
        moving= load_nifti(os.path.join(self.templatePath, self.target+ '_FA.nii.gz'))
        imgs, masks= read_caselist(self.target_csv)
        for imgPath, maskPath in (imgs, masks):

            directory= os.path.dirname(imgPath)
            inPrefix= imgPath.split('.')[0]
            prefix= os.path.split(inPrefix)[-1]
            
            b0, shm_coeff, fit_matrix= dti_harm(imgPath, maskPath, self.N_shm)

            outPrefix= os.path.join(directory, 'harm', 'ToSubjectSpace_'+ prefix)
            antsReg(imgPath, maskPath, moving, outPrefix)
            antsApply(self.templatePath, directory, prefix, self.N_shm)

            # harmonize the rish features
            mappedFile= ring_masking(directory, prefix, maskPath, self.N_shm, shm_coeff, fit_matrix, b0)
            outPrefix= os.path.join(directory, 'harm', 'rish', prefix)
            rish(mappedFile, maskPath, inPrefix, outPrefix, self.N_shm)


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


        # go through each file listed in csv, check their existence, create dti and harm
        if self.ref_csv:
            check_csv(self.ref_csv, self.force)
        check_csv(self.target_csv, self.force)



    def main(self):

        self.sanityCheck()

        if self.create:
            self.createTemplate()

        if self.process:
            self.harmonizeData()


if __name__ == '__main__':
    pipeline.run()

'''
/home/tb571/Downloads/Harmonization-Python/lib/harmonization.py \
--reference /home/tb571/Downloads/Harmonization-Python/connectom_prisma_demoData/template_debug/CaseList_MultivariateTemplate.csv \
--target /home/tb571/Downloads/Harmonization-Python/connectom_prisma_demoData/template_debug/CaseList_MultivariateTemplate1.txt \
--template /home/tb571/Downloads/Harmonization-Python/connectom_prisma_demoData/abc \
--dti /home/tb571/Downloads/Harmonization-Python/connectom_prisma_demoData/def \
--harm /home/tb571/Downloads/Harmonization-Python/connectom_prisma_demoData/ghi \
--travelHeads \
--resample \
--force \
--N_shm 6

'''