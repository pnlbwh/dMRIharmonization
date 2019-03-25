#!/usr/bin/env python

from plumbum import cli
from distutils.spawn import find_executable
import os, shutil, psutil, multiprocessing
import numpy as np

import warnings
with warnings.catch_warnings():
    warnings.filterwarnings("ignore", category=FutureWarning)
    from dipy.io.image import load_nifti

N_CPU= psutil.cpu_count()
SCRIPTDIR= os.path.dirname(__file__)

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


class pipeline(cli.Application):

    """Template creation and harmonization"""

    ref_csv = cli.SwitchAttr(
        ['--ref_list'],
        cli.ExistingFile,
        help='reference csv/txt file with first column for dwi and 2nd column for mask: dwi1,mask1\ndwi2,mask2\n...',
        mandatory=False)

    target_csv = cli.SwitchAttr(
        ['--tar_list'],
        cli.ExistingFile,
        help='target csv/txt file with first column for dwi and 2nd column for mask: dwi1,mask1\ndwi2,mask2\n...',
        mandatory=True)

    harm_csv = cli.SwitchAttr(
        ['--harm_list'],
        cli.ExistingFile,
        help='harmonized csv/txt file with first column for dwi and 2nd column for mask: dwi1,mask1\ndwi2,mask2\n...',
        mandatory=False)

    templatePath = cli.SwitchAttr(
        ['--template'],
        help='template directory',
        mandatory=True)

    N_shm = cli.SwitchAttr(
        ['--nshm'],
        help='spherical harmonic order',
        default= 6)

    N_proc = cli.SwitchAttr(
        '--nproc',
        help= 'number of processes/threads to use (-1 for all available, may slow down your system)',
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


    def createTemplate(self):

        from buildTemplate import difference_calc, antsMult, warp_bands, \
            dti_stat, rish_stat, template_masking, createAntsCaselist
        from preprocess import read_caselist, common_processing

        # check directory existence
        check_dir(self.templatePath, self.force)

        # go through each file listed in csv, check their existence, create dti and harm directories
        check_csv(self.ref_csv, self.force)
        check_csv(self.target_csv, self.force)

        # createTemplate steps -----------------------------------------------------------------------------------------

        # read image lists
        refImgs, refMasks= common_processing(self.ref_csv)
        if not self.ref_csv.endswith('.modified'):
            self.ref_csv += '.modified'
        # debug: use the following line to omit processing again
        # refImgs, refMasks = read_caselist(self.ref_csv)

        targetImgs, targetMasks= common_processing(self.target_csv)
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

        # ATTN: antsMultivariateTemplateConstruction2.sh requires '/' at the end of templatePath
        if not self.templatePath.endswith('/'):
            self.templatePath= self.templatePath+ '/'
        # ATTN: antsMultivariateTemplateConstruction2.sh requires absolute path for caselist
        antsMult(os.path.abspath(antsMultCaselist), self.templatePath)

        # load templateAffine
        templateAffine= load_nifti(os.path.join(self.templatePath, 'template0.nii.gz'))[1]


        # warp mask, dti, and rish bands
        pool = multiprocessing.Pool(self.N_proc)
        for imgPath, maskPath in zip(imgs, masks):
            pool.apply_async(func= warp_bands, args= (imgPath, maskPath, self.templatePath, ))

        pool.close()
        pool.join()

        print('dti statistics: mean, std(FA, MD) calculation of reference site')
        refMaskPath= dti_stat(self.reference, refImgs, refMasks, self.templatePath, templateAffine)
        print('dti statistics: mean, std(FA, MD) calculation of target site')
        targetMaskPath= dti_stat(self.target, targetImgs, targetMasks, self.templatePath, templateAffine)

        print('masking dti statistics of reference site')
        _= template_masking(refMaskPath, targetMaskPath, self.templatePath, self.reference)
        print('masking dti statistics of target site')
        templateMask= template_masking(refMaskPath, targetMaskPath, self.templatePath, self.target)

        print('rish_statistics mean, std(L{i}) calculation of reference site')
        rish_stat(self.reference, imgs, self.templatePath, templateAffine)
        print('rish_statistics mean, std(L{i}) calculation of target site')
        rish_stat(self.target, imgs, self.templatePath, templateAffine)

        print('calculating scale map for diffusionMeasures')
        difference_calc(self.reference, self.target, refImgs, targetImgs, self.templatePath, templateAffine,
                        'dti', templateMask, self.diffusionMeasures)

        print('calculating scale map for rishFeatures')
        difference_calc(self.reference, self.target, refImgs, targetImgs, self.templatePath, templateAffine,
                        'harm', templateMask, [f'L{i}' for i in range(0, self.N_shm+1, 2)])


        print('\n\nTemplate creation completed \n\n')


    def harmonizeData(self):

        from cleanOutliers import reconst
        from preprocess import read_caselist, dti_harm

        # check the templatePath
        if not os.path.exists(self.templatePath):
            raise NotADirectoryError(f'{self.templatePath} does not exist')
        else:
            if not os.listdir(self.templatePath):
                raise ValueError(f'{self.templatePath} is empty')

        # go through each file listed in csv, check their existence, create dti and harm directories
        check_csv(self.target_csv, self.force)

        # if self.debug:
        #     # calcuate diffusion measures of target site before any processing so we are able to compare
        #     # with the ones after harmonization
        #     pool = multiprocessing.Pool(self.N_proc)
        #     res= pool.map_async(dti_harm, zip(imgs, masks))
        #
        #     res.get()
        #     pool.close()
        #     pool.join()

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


        self.harm_csv= self.target_csv+'.harmonized'
        fh= open(self.harm_csv, 'w')
        pool = multiprocessing.Pool(self.N_proc)
        res= []
        for imgPath, maskPath in zip(imgs, masks):
            res.append(pool.apply_async(func= reconst, args= (imgPath, maskPath, moving, self.templatePath, preFlag, )))

        for r in res:
            imgPath, maskPath, harmImg, harmMask= r.get()

            if preFlag:
                fm.write(imgPath + ',' + maskPath + '\n')
            fh.write(harmImg + ',' + harmMask + '\n')


        pool.close()
        pool.join()

        if preFlag:
            fm.close()
        fh.close()
        print('\n\nHarmonization completed\n\n')


    def post_debug(self):

        from debug_fa import sub2tmp2mni, analyzeStat

        print('\n\n Reference site')
        sub2tmp2mni(self.templatePath, self.reference, self.ref_csv)
        ref_mean = analyzeStat(self.ref_csv)

        print('\n\n Target site before harmonization')
        sub2tmp2mni(self.templatePath, self.target, self.target_csv)
        target_mean_before = analyzeStat(self.target_csv)

        print('\n\n Target site after harmonization')
        sub2tmp2mni(self.templatePath, self.target, self.harm_csv)
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
        self.N_proc= int(self.N_proc)
        if self.N_proc==-1:
            self.N_proc= N_CPU

        with open(os.path.join(SCRIPTDIR,'config.ini'),'w') as f:
            f.write('[DEFAULT]\n')
            f.write(f'N_shm = {self.N_shm}\n')
            f.write(f'N_proc = {self.N_proc}\n')
            f.write(f'resample = {self.resample if self.resample else 0}\n')
            f.write(f'bvalMap = {self.bvalMap if self.bvalMap else 0}\n')
            f.write(f'denoise = {1 if self.denoise else 0}\n')
            f.write(f'travelHeads = {1 if self.travelHeads else 0}\n')
            f.write(f'debug = {1 if self.debug else 0}\n')
            f.write('diffusionMeasures = {}\n'.format((',').join(self.diffusionMeasures)))


        self.sanityCheck()

        if self.create:
            self.createTemplate()

        if self.process:
            self.harmonizeData()

        if self.debug:
            self.post_debug()


if __name__ == '__main__':
    pipeline.run()

