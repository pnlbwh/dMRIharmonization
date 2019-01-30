#!/usr/bin/env python

import warnings, os
with warnings.catch_warnings():
    warnings.filterwarnings("ignore", category=FutureWarning)
    import nibabel as nib
    from dipy.segment.mask import applymask

from plumbum.cmd import antsApplyTransforms
from plumbum import FG
from buildTemplate import read_caselist
from cleanOutliers import antsReg

module_dir= os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..',))

def sub2tmp2mni(templatePath, siteName, caselist, diffusionMeasures):

    # obtain the transform
    moving = os.path.join(templatePath, f'Mean_{siteName}_FA.nii.gz')
    fixed = os.path.join(module_dir, 'IITAtlas', 'IITmean_FA.nii.gz')

    outPrefix= os.path.join(templatePath, f'TemplateToMNI_{siteName}')
    warp2mni= outPrefix+'1Warp.nii.gz'
    trans2mni= outPrefix+'0GenericAffine.mat'
    # if not os.path.exists(warp2mni):
    antsReg(fixed, None, moving, outPrefix)

    imgs, _= read_caselist(caselist)

    # TODO: parellelize
    ref = False
    tar_harm = False
    tar_unharm = False
    for imgPath in imgs:

        print(f'Warping {imgPath} diffusion measures to standard space')
        directory= os.path.dirname(imgPath)
        inPrefix= imgPath.split('.')[0]
        prefix= os.path.split(inPrefix)[-1]

        # reference sites have been already warped to template space in buildTemplate.py: warp_bands()
        # warped data are os.path.join(directory, 'dti', prefix + f'_WarpedFA.nii.gz')

        # template was warped to target sites in cleanOutliers.py: antsApply()
        # to go from target sites to template space, use the following files:
        trans2sub= os.path.join(directory, 'harm', f'ToSubjectSpace_{prefix}0GenericAffine.mat')
        warp2tmp= os.path.join(directory, 'harm', f'ToSubjectSpace_{prefix}1InverseWarp.nii.gz')

        # only harmonized data needs to SyN register to template space, identified with 'harmonized_' in prefix
        if 'harmonized_' in prefix:
            dmImg = os.path.join(directory, 'dti', prefix + f'_FA.nii.gz')
            dmTmp = os.path.join(templatePath, f'Mean_{siteName}_FA.nii.gz')
            maskTmp = os.path.join(templatePath, f'{siteName}_Mask.nii.gz')
            outPrefix = os.path.join(directory, 'dti', prefix+ '_FA')
            # antsReg(img, mask, mov, outPrefix)
            antsReg(dmTmp, maskTmp, dmImg, outPrefix)
            warp2tmp = outPrefix + '1Warp.nii.gz'
            trans2tmp = outPrefix + '0GenericAffine.mat'
            tar_harm= True

        elif os.path.exists(trans2sub) and os.path.exists(warp2tmp):
            tar_unharm= True

        else:
            ref= True


        for dm in diffusionMeasures:

            output= os.path.join(directory, 'dti', prefix+ f'_InMNI_{dm}.nii.gz')

            if ref:
                moving = os.path.join(directory, 'dti', prefix + f'_Warped{dm}.nii.gz')
                # warp diffusion measure to MNI space directly
                antsApplyTransforms[
                    '-d', '3',
                    '-i', moving,
                    '-o', output,
                    '-r', fixed,
                    '-t', warp2mni, trans2mni
                ] & FG

            elif tar_unharm:
                moving = os.path.join(directory, 'dti', prefix + f'_{dm}.nii.gz')
                # warp diffusion measure to template space first, then to MNI space
                antsApplyTransforms[
                    '-d', '3',
                    '-i', moving,
                    '-o', output,
                    '-r', fixed,
                    '-t', warp2mni, trans2mni, f'[{trans2sub},1]', warp2tmp
                    ] & FG

            elif tar_harm:
                moving = os.path.join(directory, 'dti', prefix + f'_{dm}.nii.gz')
                # warp diffusion measure to template space first, then to MNI space
                antsApplyTransforms[
                    '-d', '3',
                    '-i', moving,
                    '-o', output,
                    '-r', fixed,
                    '-t', warp2mni, trans2mni, warp2tmp, trans2tmp
                    ] & FG


def analyzeStat(file):
    '''
    :param file: list of (FA or MD or GFA) that are already in MNI space
    :return: mean of the images
    '''

    skel= nib.load(os.path.join(module_dir, 'IITAtlas', 'IITmean_FA_skeleton.nii.gz'))
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