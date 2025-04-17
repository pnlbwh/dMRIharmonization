#!/usr/bin/env python

# ===============================================================================
# dMRIharmonization (2018) pipeline is written by-
#
# TASHRIF BILLAH
# Brigham and Women's Hospital/Harvard Medical School
# tbillah@bwh.harvard.edu, tashrifbillah@gmail.com
#
# ===============================================================================
# See details at https://github.com/pnlbwh/dMRIharmonization
# Submit issues at https://github.com/pnlbwh/dMRIharmonization/issues
# View LICENSE at https://github.com/pnlbwh/dMRIharmonization/blob/master/LICENSE
# ===============================================================================

from plumbum import cli
from shutil import which
import multiprocessing
import io

from determineNshm import verifyNshmForAll, determineNshm
from util import *
from fileUtil import read_caselist, check_dir, check_csv

N_CPU= multiprocessing.cpu_count()
SCRIPTDIR= dirname(__file__)


class pipeline(cli.Application):

    """
    ===============================================================================
    dMRIharmonization (2018) pipeline is written by-

    TASHRIF BILLAH
    Brigham and Women's Hospital/Harvard Medical School
    tbillah@bwh.harvard.edu, tashrifbillah@gmail.com

    ===============================================================================
    See details at https://github.com/pnlbwh/dMRIharmonization
    Submit issues at https://github.com/pnlbwh/dMRIharmonization/issues
    View LICENSE at https://github.com/pnlbwh/dMRIharmonization/blob/master/LICENSE
    ===============================================================================

    Template creation, harmonization, and debugging
    """

    ref_csv = cli.SwitchAttr(
        ['--ref_list'],
        cli.ExistingFile,
        help='reference csv/txt file with first column for dwi and 2nd column for mask: dwi1,mask1\\n dwi2,mask2\\n...',
        mandatory=False)

    target_csv = cli.SwitchAttr(
        ['--tar_list'],
        cli.ExistingFile,
        help='target csv/txt file with first column for dwi and 2nd column for mask: dwi1,mask1\\n dwi2,mask2\\n...',
        mandatory=False)

    harm_csv = cli.SwitchAttr(
        ['--harm_list'],
        cli.ExistingFile,
        help='harmonized csv/txt file with first column for dwi and 2nd column for mask: dwi1,mask1\\n dwi2,mask2\\n...',
        mandatory=False)

    templatePath = cli.SwitchAttr(
        ['--template'],
        help='template directory',
        mandatory=True)

    N_shm = cli.SwitchAttr(
        ['--nshm'],
        help='spherical harmonic order, by default maximum possible is used',
        default= '-1')

    N_proc = cli.SwitchAttr(
        '--nproc',
        help= 'number of processes/threads to use (-1 for all available, may slow down your system)',
        default= 4)

    N_zero = cli.SwitchAttr(
        '--nzero',
        help= 'number of zero padding for denoising skull region during signal reconstruction',
        default= 10)

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

    bshell_b = cli.SwitchAttr(
        '--bshell_b',
        help='bvalue of the bshell, needed for multi-shell data only',
        mandatory= False)

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
        '--ref_name',
        help= 'reference site name',
        mandatory= False)

    target = cli.SwitchAttr(
        '--tar_name',
        help= 'target site name',
        mandatory= True)

    verbose= cli.Flag(
        '--verbose',
        help='print everything to STDOUT',
        default= False)


    diffusionMeasures = ['MD', 'FA', 'GFA']


    def createTemplate(self):

        from buildTemplate import difference_calc, antsMult, warp_bands, \
            dti_stat, rish_stat, template_masking, createAntsCaselist
        from preprocess import common_processing

        # check directory existence
        check_dir(self.templatePath, self.force)

        # go through each file listed in csv, check their existence, create dti and harm directories
        check_csv(self.ref_csv, self.force)
        check_csv(self.target_csv, self.force)

        # createTemplate steps -----------------------------------------------------------------------------------------

        # read image lists
        refImgs, refMasks= common_processing(self.ref_unproc_csv)
        if not self.ref_csv.endswith('.modified'):
            self.ref_csv += '.modified'
        # debug: use the following line to omit processing again
        # refImgs, refMasks = read_caselist(self.ref_csv)

        targetImgs, targetMasks= common_processing(self.tar_unproc_csv)
        if not self.target_csv.endswith('.modified'):
            self.target_csv += '.modified'
        # debug: use the following line to omit processing again
        # targetImgs, targetMasks = read_caselist(self.target_csv)

        imgs= refImgs+targetImgs
        masks= refMasks+targetMasks

        # create caselist for antsMult
        antsMultCaselist= pjoin(self.templatePath, 'antsMultCaselist.txt')
        createAntsCaselist(imgs, antsMultCaselist)

        # run ANTS multivariate template construction

        # ATTN: antsMultivariateTemplateConstruction2.sh requires '/' at the end of templatePath
        if not self.templatePath.endswith('/'):
            self.templatePath += '/'
        
        # check if template was created earlier
        prevTemplateFile= pjoin(self.templatePath, 'prevTemplateCompletion')
        template0= pjoin(self.templatePath, 'template0.nii.gz')
        if not isfile(prevTemplateFile):
            # ATTN: antsMultivariateTemplateConstruction2.sh requires absolute path for caselist
            antsMult(abspath(antsMultCaselist), self.templatePath)
        else:
            warnings.warn(f'Using {template0} which was created before')
            
        # load templateHdr
        templateHdr= load(template0).header


        # warp mask, dti, and rish bands
        if self.N_proc==1:
            for imgPath, maskPath in zip(imgs, masks):
                warp_bands(imgPath, maskPath, self.templatePath)
        
        elif self.N_proc>1:
            pool = multiprocessing.Pool(self.N_proc)
            for imgPath, maskPath in zip(imgs, masks):
                pool.apply_async(func= warp_bands, args= (imgPath, maskPath, self.templatePath,))

            pool.close()
            pool.join()

        print('calculating dti statistics i.e. mean, std for reference site')
        refMaskPath= dti_stat(self.reference, refImgs, refMasks, self.templatePath, templateHdr)
        print('calculating dti statistics i.e. mean, std for target site')
        targetMaskPath= dti_stat(self.target, targetImgs, targetMasks, self.templatePath, templateHdr)

        print('masking dti statistics of reference site')
        _= template_masking(refMaskPath, targetMaskPath, self.templatePath, self.reference)
        print('masking dti statistics of target site')
        templateMask= template_masking(refMaskPath, targetMaskPath, self.templatePath, self.target)

        print('calculating rish_statistics i.e. mean, std calculation of reference site')
        rish_stat(self.reference, refImgs, self.templatePath, templateHdr)
        print('calculating rish_statistics i.e. mean, std calculation of target site')
        rish_stat(self.target, targetImgs, self.templatePath, templateHdr)

        print('calculating templates for diffusionMeasures')
        difference_calc(self.reference, self.target, refImgs, targetImgs, self.templatePath, templateHdr,
                        templateMask, self.diffusionMeasures)

        print('calculating templates for rishFeatures')
        difference_calc(self.reference, self.target, refImgs, targetImgs, self.templatePath, templateHdr,
                        templateMask, [f'L{i}' for i in range(0, self.N_shm+1, 2)])

        
        # write a flag in templatePath that can be used to see if template was created earlier
        if not isfile(prevTemplateFile):
            with open(prevTemplateFile, 'w'):
                pass

        print('\n\nTemplate creation completed \n\n')


    def harmonizeData(self):

        from reconstSignal import reconst, approx
        from preprocess import dti_harm, common_processing, preprocessing

        # check the templatePath
        if not exists(self.templatePath):
            raise NotADirectoryError(f'{self.templatePath} does not exist')
        else:
            if not listdir(self.templatePath):
                raise ValueError(f'{self.templatePath} is empty')



        # fit spherical harmonics on reference site
        if self.debug and self.ref_csv:
            check_csv(self.ref_unproc_csv, self.force)
            refImgs, refMasks= read_caselist(self.ref_unproc_csv)

            if self.N_proc==1:
                attributes=[]
                for imgPath, maskPath in zip(refImgs, refMasks):
                    attributes.append(preprocessing(imgPath, maskPath))
            
            elif self.N_proc>1:
                res= []
                pool = multiprocessing.Pool(self.N_proc)
                for imgPath, maskPath in zip(refImgs, refMasks):
                    res.append(pool.apply_async(func=preprocessing, args=(imgPath, maskPath)))

                attributes = [r.get() for r in res]

                pool.close()
                pool.join()

            for i in range(len(refImgs)):
                refImgs[i] = attributes[i][0]
                refMasks[i] = attributes[i][1]

            if self.N_proc==1:
                for imgPath, maskPath in zip(refImgs, refMasks):
                    approx(imgPath,maskPath)

            elif self.N_proc>1:
                pool = multiprocessing.Pool(self.N_proc)
                for imgPath, maskPath in zip(refImgs, refMasks):
                    pool.apply_async(func= approx, args=(imgPath,maskPath,))

                pool.close()
                pool.join()



        # go through each file listed in csv, check their existence, create dti and harm directories
        check_csv(self.target_csv, self.force)
        targetImgs, targetMasks= common_processing(self.tar_unproc_csv)


        # reconstSignal steps ------------------------------------------------------------------------------------------

        # read target image list
        moving= pjoin(self.templatePath, f'Mean_{self.target}_FA.nii.gz')

        if not self.target_csv.endswith('.modified'):
            self.target_csv += '.modified'


        self.harm_csv= self.target_csv+'.harmonized'
        fh= open(self.harm_csv, 'w')
        if self.N_proc==1:
            res=[]
            for imgPath, maskPath in zip(targetImgs, targetMasks):
                res.append(reconst(imgPath, maskPath, moving, self.templatePath))

            for r in res:
                harmImg, harmMask= r
                fh.write(harmImg + ',' + harmMask + '\n')

        elif self.N_proc>1:
            pool = multiprocessing.Pool(self.N_proc)
            res= []
            for imgPath, maskPath in zip(targetImgs, targetMasks):
                res.append(pool.apply_async(func= reconst, args= (imgPath, maskPath, moving, self.templatePath,)))

            for r in res:
                harmImg, harmMask= r.get()
                fh.write(harmImg + ',' + harmMask + '\n')


            pool.close()
            pool.join()

        fh.close()
        
        
        if self.debug:
            harmImgs, harmMasks= read_caselist(self.harm_csv)

            if self.N_proc==1:
                for imgPath,maskPath in zip(harmImgs,harmMasks):
                    dti_harm(imgPath,maskPath)

            elif self.N_proc>=1:
                pool = multiprocessing.Pool(self.N_proc)
                for imgPath,maskPath in zip(harmImgs,harmMasks):
                    pool.apply_async(func= dti_harm, args= (imgPath,maskPath,))
                pool.close()
                pool.join()
            
        print('\n\nHarmonization completed\n\n')


    def post_debug(self):

        from debug_fa import sub2tmp2mni

        print('\n\n Reference site')
        sub2tmp2mni(self.templatePath, self.reference, self.ref_csv, ref= True)

        print('\n\n Target site before harmonization')
        sub2tmp2mni(self.templatePath, self.target, self.tar_unproc_csv, tar_unproc= True)

        print('\n\n Target site after harmonization')
        sub2tmp2mni(self.templatePath, self.target, self.harm_csv, tar_harm= True)


        self.showStat()


    def showStat(self):

        from debug_fa import analyzeStat
        from datetime import datetime
        from harm_plot import generate_csv, harm_plot
        import pandas as pd
        
        print('\n\nComputing statistics:')
        
        print(f'{self.reference} site')
        ref_mean = analyzeStat(self.ref_csv, self.templatePath)
        generate_csv(self.ref_csv, ref_mean, pjoin(self.templatePath, self.reference))

        print(f'{self.target} site before harmonization')
        target_mean_before = analyzeStat(self.tar_unproc_csv, self.templatePath)
        generate_csv(self.tar_unproc_csv, target_mean_before, pjoin(self.templatePath, self.target+'_before'))

        print(f'{self.target} site after harmonization')
        target_mean_after = analyzeStat(self.harm_csv, self.templatePath)
        generate_csv(self.harm_csv, target_mean_after, pjoin(self.templatePath, self.target+'_after'))

        
        print('\n\nPrinting statistics:')
        # save statistics for future
        statFile= pjoin(self.templatePath, 'meanFAstat.csv')
        if isfile(statFile):
            df= pd.read_csv(statFile)
        else:
            timestamp= datetime.now().strftime('%m/%d/%y %H:%M')
            sites= [f'{self.reference}',f'{self.target}_before',f'{self.target}_after']
            df= pd.DataFrame({timestamp:sites})

        header= f'mean meanFA b{self.bshell_b}'
        value= [np.mean(x) for x in [ref_mean, target_mean_before, target_mean_after]]
        df= df.assign(**{header:value})
        df.to_csv(statFile, index=False)
        

        # print statistics on console
        with open(statFile) as f:
            print(f.read())
            
        
        # generate graph
        ebar= harm_plot([ref_mean, target_mean_before, target_mean_after],
                         labels=[self.reference, self.target+'_before', self.target+'_after'],
                         outPrefix=pjoin(self.templatePath,'meanFAstat'))

        print(f'\nDetailed statistics, summary results, and demonstrative plots are saved in:\n\n{self.templatePath}/*_stat.csv'
              f'\n{statFile}\n{ebar}\n')


    def sanityCheck(self):

        if not (self.create or self.process or self.debug):
            raise AttributeError('No option selected, ' 
                                 'specify one (or many of) creation, harmonization, and debug flags')

        # check ants commands
        external_commands= [
            'antsMultivariateTemplateConstruction2.sh',
            'antsApplyTransforms',
            'antsRegistrationSyNQuick.sh',
            'unring.a64']

        for cmd in external_commands:
            exe= which(cmd)
            if not exe:
                raise EnvironmentError(f'{cmd} not found')



    def main(self):
    
        self.sanityCheck()

        self.templatePath= abspath(self.templatePath)
        self.N_shm= int(self.N_shm)
        self.N_proc= int(self.N_proc)
        if self.N_proc==-1:
            self.N_proc= N_CPU

        if self.ref_csv:
            self.ref_unproc_csv= self.ref_csv.strip('.modified')
        self.tar_unproc_csv= self.target_csv.strip('.modified')


        # check appropriateness of N_shm
        if self.N_shm!=-1 and (self.N_shm<2 or self.N_shm>8):
            raise ValueError('2<= --nshm <=8')


        # determine N_shm in default mode during template creation
        if self.N_shm==-1 and self.create:
            if self.ref_csv:
                ref_nshm_img = read_caselist(self.ref_csv)[0][0]
            elif self.target_csv:
                ref_nshm_img = read_caselist(self.target_csv)[0][0]

            directory= dirname(ref_nshm_img)
            prefix= basename(ref_nshm_img).split('.nii')[0]
            bvalFile= pjoin(directory, prefix+'.bval')
            self.N_shm, _= determineNshm(bvalFile)


        # automatic determination of N_shm during data harmonization is limited by N_shm used during template creation
        # Scale_L{i}.nii.gz of <= {N_shm during template creation} are present only
        elif self.N_shm==-1 and self.process:
            for i in range(0,8,2):
                if isfile(pjoin(self.templatePath, f'Scale_L{i}.nii.gz')):
                    self.N_shm= i
                else:
                    break


        # verify validity of provided/determined N_shm for all subjects
        # single-shell-ness is verified inside verifyNshmForAll
        if self.ref_csv:
            verifyNshmForAll(self.ref_csv, self.N_shm)
        if self.target_csv:
            verifyNshmForAll(self.target_csv, self.N_shm)


        # write config file to temporary directory
        configFile= pjoin(gettempdir(),f'harm_config_{getpid()}.ini')
        with open(configFile,'w') as f:
            f.write('[DEFAULT]\n')
            f.write(f'N_shm = {self.N_shm}\n')
            f.write(f'N_proc = {self.N_proc}\n')
            f.write(f'N_zero = {self.N_zero}\n')
            f.write(f'resample = {self.resample if self.resample else 0}\n')
            f.write(f'bvalMap = {self.bvalMap if self.bvalMap else 0}\n')
            f.write(f'bshell_b = {self.bshell_b}\n')
            f.write(f'denoise = {1 if self.denoise else 0}\n')
            f.write(f'travelHeads = {1 if self.travelHeads else 0}\n')
            f.write(f'debug = {1 if self.debug else 0}\n')
            f.write(f'force = {1 if self.force else 0}\n')
            f.write(f'verbose = {1 if self.verbose else 0}\n')
            f.write('diffusionMeasures = {}\n'.format((',').join(self.diffusionMeasures)))


        if self.create:
            self.createTemplate()
            import fileinput
            for line in fileinput.input(configFile, inplace=True):
                if 'force' in line:
                    print('force = 0')
                else:
                    print(line)
            self.force= False
            
        if self.process:
            self.harmonizeData()

        if self.create and self.process and self.debug:
            self.post_debug()


        remove(configFile)


if __name__ == '__main__':
    pipeline.run()

