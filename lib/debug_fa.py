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

from plumbum.cmd import antsApplyTransforms
from plumbum import FG
from preprocess import read_caselist
from reconstSignal import antsReg
import multiprocessing
from util import *

SCRIPTDIR= dirname(__file__)
ROOTDIR= abspath(pjoin(SCRIPTDIR, '..'))
mniTmp = pjoin(ROOTDIR, 'IITAtlas', 'IITmean_FA.nii.gz')

config = ConfigParser()
config.read(pjoin(gettempdir(),f'harm_config_{getpid()}.ini'))
N_proc = int(config['DEFAULT']['N_proc'])
diffusionMeasures = [x for x in config['DEFAULT']['diffusionMeasures'].split(',')]

def register_reference(imgPath, warp2mni, trans2mni, templatePath):

    print(f'Warping {imgPath} diffusion measures to standard space')
    directory = dirname(imgPath)
    inPrefix = imgPath.split('.nii')[0]
    prefix = psplit(inPrefix)[-1]

    for dm in diffusionMeasures:

        output = pjoin(templatePath, prefix + f'_InMNI_{dm}.nii.gz')

        # reference site have been already warped to reference template space in buildTemplate.py: warp_bands()
        # warped data are pjoin(directory, 'dti', prefix + f'_WarpedFA.nii.gz')
        moving = pjoin(templatePath, prefix + f'_Warped{dm}.nii.gz')

        # so warp diffusion measures that are in template space to MNI space
        antsApplyTransforms[
            '-d', '3',
            '-i', moving,
            '-o', output,
            '-r', mniTmp,
            '-t', warp2mni, trans2mni
        ] & FG


def register_target(imgPath, templatePath):

    print(f'Warping {imgPath} diffusion measures to standard space')
    directory = dirname(imgPath)
    inPrefix = imgPath.split('.nii')[0]
    prefix = psplit(inPrefix)[-1]

    dmImg = pjoin(directory, 'dti', prefix + f'_FA.nii.gz')

    outPrefix = pjoin(templatePath, prefix + '_FA_ToMNI_')
    warp2mni = outPrefix + '1Warp.nii.gz'
    trans2mni = outPrefix + '0GenericAffine.mat'
    # unprocessed target data is given, so in case multiple debug is needed, pass the registration
    if not exists(warp2mni):
        antsReg(mniTmp, None, dmImg, outPrefix)

    for dm in diffusionMeasures:
        output = pjoin(templatePath, prefix + f'_InMNI_{dm}.nii.gz')

        moving = pjoin(directory, 'dti', prefix + f'_{dm}.nii.gz')
        # warp diffusion measures to MNI space directly
        antsApplyTransforms[
            '-d', '3',
            '-i', moving,
            '-o', output,
            '-r', mniTmp,
            '-t', warp2mni, trans2mni,
        ] & FG


def register_harmonized(imgPath, warp2mni, trans2mni, templatePath, siteName):

    print(f'Warping {imgPath} diffusion measures to standard space')
    directory = dirname(imgPath)
    inPrefix = imgPath.split('.nii')[0]
    prefix = psplit(inPrefix)[-1]

    dmImg = pjoin(directory, 'dti', prefix + f'_FA.nii.gz')
    dmTmp = pjoin(templatePath, f'Mean_{siteName}_FA.nii.gz')
    maskTmp = pjoin(templatePath, f'{siteName}_Mask.nii.gz')
    outPrefix = pjoin(templatePath, prefix + '_FA_ToMNI')
    warp2tmp = outPrefix + '1Warp.nii.gz'
    trans2tmp = outPrefix + '0GenericAffine.mat'
    # signal reconstruction might change with zero padding size, median filtering kernel size, and harmonized mask
    # so in case multiple debug is needed, redo the registration
    antsReg(dmTmp, maskTmp, dmImg, outPrefix)

    for dm in diffusionMeasures:
        output = pjoin(templatePath, prefix + f'_{dm}_ToTmpWarped.nii.gz')

        moving = pjoin(directory, 'dti', prefix + f'_{dm}.nii.gz')
        # warp diffusion measures to template space first, then to MNI space
        antsApplyTransforms[
            '-d', '3',
            '-i', moving,
            '-o', output,
            '-r', dmTmp,
            '-t', warp2tmp, trans2tmp
        ] & FG

        output = pjoin(templatePath, prefix + f'_InMNI_{dm}.nii.gz')

        moving = pjoin(templatePath, prefix + f'_{dm}_ToTmpWarped.nii.gz')


        antsApplyTransforms[
            '-d', '3',
            '-i', moving,
            '-o', output,
            '-r', mniTmp,
            '-t', warp2mni, trans2mni
        ] & FG


def sub2tmp2mni(templatePath, siteName, caselist, ref= False, tar_unproc= False, tar_harm= False):

    # obtain the transform
    moving = pjoin(templatePath, f'Mean_{siteName}_FA.nii.gz')

    outPrefix= pjoin(templatePath, f'TemplateToMNI_{siteName}')
    warp2mni= outPrefix+'1Warp.nii.gz'
    trans2mni= outPrefix+'0GenericAffine.mat'
    # template is created once, it is expected that the user wants to keep the template same during debugging
    # so in case multiple debug is needed, pass the registration
    if not exists(warp2mni):
        antsReg(mniTmp, None, moving, outPrefix)

    imgs, _= read_caselist(caselist)

    pool= multiprocessing.Pool(N_proc)
    for imgPath in imgs:

        if ref:
            pool.apply_async(func= register_reference, args= (imgPath, warp2mni, trans2mni, templatePath, ))
        elif tar_unproc:
            pool.apply_async(func= register_target, args= (imgPath, templatePath, ))
        elif tar_harm:
            pool.apply_async(func= register_harmonized, args= (imgPath, warp2mni, trans2mni, templatePath, siteName, ))

    pool.close()
    pool.join()


def analyzeStat(file, templatePath):

    skel= load(pjoin(ROOTDIR, 'IITAtlas', 'IITmean_FA_skeleton.nii.gz'))
    skel_mask= (skel.get_data()>0)*1.

    imgs, _ = read_caselist(file)

    meanAttr=[]
    for imgPath in imgs:
        inPrefix = imgPath.split('.nii')[0]
        prefix = psplit(inPrefix)[-1]

        faImg= pjoin(templatePath, prefix + f'_InMNI_FA.nii.gz')
        data= load(faImg).get_data()
        temp= data*skel_mask
        meanAttr.append(temp[temp>0].mean())

    return meanAttr
