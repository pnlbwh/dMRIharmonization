#!/usr/bin/env python

import warnings, os, configparser, multiprocessing
with warnings.catch_warnings():
    warnings.filterwarnings("ignore", category=FutureWarning)
    import nibabel as nib

from plumbum.cmd import antsApplyTransforms
from plumbum import FG
from preprocess import read_caselist
from cleanOutliers import antsReg

SCRIPTDIR= os.path.dirname(__file__)
ROOTDIR= os.path.abspath(os.path.join(SCRIPTDIR, '..'))
mniTmp = os.path.join(ROOTDIR, 'IITAtlas', 'IITmean_FA.nii.gz')

config = configparser.ConfigParser()
config.read(os.path.join(SCRIPTDIR,'config.ini'))
N_proc = int(config['DEFAULT']['N_proc'])
diffusionMeasures = [x for x in config['DEFAULT']['diffusionMeasures'].split(',')]

def register_reference(imgPath, mniTmp, warp2mni, trans2mni):

    print(f'Warping {imgPath} diffusion measures to standard space')
    directory = os.path.dirname(imgPath)
    inPrefix = imgPath.split('.')[0]
    prefix = os.path.split(inPrefix)[-1]

    for dm in diffusionMeasures:

        output = os.path.join(directory, 'dti', prefix + f'_InMNI_{dm}.nii.gz')

        # reference site have been already warped to reference template space in buildTemplate.py: warp_bands()
        # warped data are os.path.join(directory, 'dti', prefix + f'_WarpedFA.nii.gz')
        moving = os.path.join(directory, 'dti', prefix + f'_Warped{dm}.nii.gz')

        # so warp diffusion measure to MNI space directly
        antsApplyTransforms[
            '-d', '3',
            '-i', moving,
            '-o', output,
            '-r', mniTmp,
            '-t', warp2mni, trans2mni
        ] & FG


def register_target(imgPath, mniTmp):

    print(f'Warping {imgPath} diffusion measures to standard space')
    directory = os.path.dirname(imgPath)
    inPrefix = imgPath.split('.')[0]
    prefix = os.path.split(inPrefix)[-1]

    dmImg = os.path.join(directory, 'dti', prefix + f'_FA.nii.gz')

    outPrefix = os.path.join(directory, 'dti', prefix + '_FA_ToMNI_')
    warp2mni = outPrefix + '1Warp.nii.gz'
    trans2mni = outPrefix + '0GenericAffine.mat'
    # unprocessed target data is given, so in case multiple debug is needed, pass the registration
    if not os.path.exists(warp2mni):
        antsReg(mniTmp, None, dmImg, outPrefix)

    for dm in diffusionMeasures:
        output = os.path.join(directory, 'dti', prefix + f'_InMNI_{dm}.nii.gz')

        moving = os.path.join(directory, 'dti', prefix + f'_{dm}.nii.gz')
        # warp diffusion measure to template space first, then to MNI space
        antsApplyTransforms[
            '-d', '3',
            '-i', moving,
            '-o', output,
            '-r', mniTmp,
            '-t', warp2mni, trans2mni,
        ] & FG


def register_harmonized(imgPath, mniTmp, warp2mni, trans2mni, templatePath, siteName):

    print(f'Warping {imgPath} diffusion measures to standard space')
    directory = os.path.dirname(imgPath)
    inPrefix = imgPath.split('.')[0]
    prefix = os.path.split(inPrefix)[-1]

    dmImg = os.path.join(directory, 'dti', prefix + f'_FA.nii.gz')
    dmTmp = os.path.join(templatePath, f'Mean_{siteName}_FA.nii.gz')
    maskTmp = os.path.join(templatePath, f'{siteName}_Mask.nii.gz')
    outPrefix = os.path.join(directory, 'dti', prefix + '_FA')
    warp2tmp = outPrefix + '1Warp.nii.gz'
    trans2tmp = outPrefix + '0GenericAffine.mat'
    # signal reconstruction might change with zero padding size, median filtering kernel size, and harmonized mask
    # so in case multiple debug is needed, redo the registration
    antsReg(dmTmp, maskTmp, dmImg, outPrefix)

    for dm in diffusionMeasures:
        output = os.path.join(directory, 'dti', prefix + f'_InMNI_{dm}.nii.gz')

        moving = os.path.join(directory, 'dti', prefix + f'_{dm}.nii.gz')
        # warp diffusion measure to template space first, then to MNI space
        antsApplyTransforms[
            '-d', '3',
            '-i', moving,
            '-o', output,
            '-r', mniTmp,
            '-t', warp2mni, trans2mni, warp2tmp, trans2tmp
        ] & FG


def sub2tmp2mni(templatePath, siteName, caselist, ref= False, tar_unproc= False, tar_harm= False):

    # obtain the transform
    moving = os.path.join(templatePath, f'Mean_{siteName}_FA.nii.gz')

    outPrefix= os.path.join(templatePath, f'TemplateToMNI_{siteName}')
    warp2mni= outPrefix+'1Warp.nii.gz'
    trans2mni= outPrefix+'0GenericAffine.mat'
    # template is created once, it is expected that the user wants to keep the template same during debugging
    # so in case multiple debug is needed, pass the registration
    if not os.path.exists(warp2mni):
        antsReg(mniTmp, None, moving, outPrefix)

    imgs, _= read_caselist(caselist)

    pool= multiprocessing.Pool(N_proc)
    for imgPath in imgs:

        if ref:
            pool.apply_async(func= register_reference, args= (imgPath, mniTmp, warp2mni, trans2mni, ))
        elif tar_unproc:
            pool.apply_async(func= register_target, args= (imgPath, mniTmp, ))
        elif tar_harm:
            pool.apply_async(func= register_harmonized, args= (imgPath, warp2mni, trans2mni, mniTmp, templatePath, siteName, ))

    pool.close()
    pool.join()


def analyzeStat(file):
    '''
    :param file: list of (FA or MD or GFA) that are already in MNI space
    :return: mean of the images
    '''

    skel= nib.load(os.path.join(ROOTDIR, 'IITAtlas', 'IITmean_FA_skeleton.nii.gz'))
    skel_mask= (skel.get_data()>0)*1.

    imgs, _ = read_caselist(file)

    meanAttr=[]
    for imgPath in imgs:
        directory = os.path.dirname(imgPath)
        inPrefix = imgPath.split('.')[0]
        prefix = os.path.split(inPrefix)[-1]

        faImg= os.path.join(directory, 'dti', prefix + f'_InMNI_FA.nii.gz')
        data= nib.load(faImg).get_data()
        temp= data*skel_mask
        meanAttr.append(temp[temp>0].mean())

    return meanAttr
