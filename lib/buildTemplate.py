#!/usr/bin/env python

import numpy as np
from subprocess import check_call
from plumbum.cmd import antsApplyTransforms
from plumbum import FG
import psutil, os
N_CPU= psutil.cpu_count()
from glob import glob
import warnings
with warnings.catch_warnings():
    warnings.filterwarnings("ignore", category=FutureWarning)
    from dipy.io.image import load_nifti, save_nifti
    from dipy.segment.mask import applymask

from scipy.ndimage import binary_opening, generate_binary_structure
from scipy.ndimage.filters import gaussian_filter
eps= 2.2204e-16


def modifiedFile(file, subDir, ext):
    prefix= os.path.basename(file).split('.')[0]
    return os.path.join(os.path.dirname(file), subDir, prefix+ext)


def applyXform(inImg, refImg, warp, trans, outImg):

    antsApplyTransforms[
        '-d', '3',
        '-i', inImg,
        '-o', outImg,
        '-r', refImg,
        '-t', warp, trans
        ] & FG


# def warp_bands(dtiPath, rishPath, maskPath, templatePath, N_shm, diffusionMeasures):
def warp_bands(imgPath, maskPath, templatePath, N_shm, diffusionMeasures):

    prefix= os.path.basename(imgPath).split('.')[0]
    directory= os.path.dirname(imgPath)
    warp = glob(os.path.join(templatePath, prefix + f'_FA*[!Inverse]Warp.nii.gz'))
    trans = glob(os.path.join(templatePath, prefix + f'_FA*GenericAffine.mat'))

    # warping the mask
    applyXform(maskPath,
               os.path.join(templatePath, 'template0.nii.gz'),
               warp, trans,
               maskPath.split('.')[0]+ 'Warped.nii.gz')


    # warping the rish features
    for i in range(0, N_shm+1, 2):
        applyXform(os.path.join(directory, 'harm', f'{prefix}_L{i}.nii.gz'),
           os.path.join(templatePath, 'template0.nii.gz'),
           warp, trans,
           os.path.join(directory, 'harm', f'{prefix}_WarpedL{i}.nii.gz'))


    # warping the diffusion measures
    for dm in diffusionMeasures:
        applyXform(os.path.join(directory, 'dti', f'{prefix}_{dm}.nii.gz'),
                   os.path.join(templatePath, 'template0.nii.gz'),
                   warp, trans,
                   os.path.join(directory, 'dti', f'{prefix}_Warped{dm}.nii.gz'))


def createAntsCaselist(imgs, file):

    with open(file,'w') as f:
        for imgPath in imgs:
            prefix= os.path.basename(imgPath).split('.')[0]
            directory= os.path.dirname(imgPath)

            FA= os.path.join(directory,'dti', f'{prefix}_FA.nii.gz')
            L0= os.path.join(directory,'harm', f'{prefix}_L0.nii.gz')
            f.write(f'{FA},{L0}\n')


def antsMult(caselist, outPrefix):

    check_call((' ').join(['antsMultivariateTemplateConstruction2.sh',
                           '-d', '3',
                           '-g', '0.2',
                           '-k', '2',
                           '-t', "BSplineSyN[0.1,26,0]",
                           '-r', '1',
                           '-c', '2',
                           '-j', str(N_CPU),
                           '-f', '8x4x2x1',
                           '-o', outPrefix,
                           caselist]), shell= True)



def dti_stat(siteName, imgs, masks, templatePath, templateAffine, diffusionMeasures):

    maskData = []
    for maskPath in masks:
        maskData.append(load_nifti(maskPath.split('.')[0]+ 'Warped.nii.gz')[0])

    morphed_mask= binary_opening(np.mean(maskData, axis= 0)>0.5, structure= generate_binary_structure(3,1))*1
    morphed_mask_name= os.path.join(templatePath, f'{siteName}_Mask.nii.gz')
    save_nifti(morphed_mask_name, morphed_mask, templateAffine)

    imgData= []
    for dm in diffusionMeasures:
        for imgPath in imgs:
            prefix = os.path.basename(imgPath).split('.')[0]
            directory = os.path.dirname(imgPath)
            imgData.append(load_nifti(os.path.join(directory, 'dti', f'{prefix}_Warped{dm}.nii.gz'))[0])

        save_nifti(os.path.join(templatePath, f'Mean_{siteName}_{dm}.nii.gz'),
                                np.mean(imgData, axis= 0), templateAffine)

        save_nifti(os.path.join(templatePath, f'Std_{siteName}_{dm}.nii.gz'),
                                np.std(imgData, axis= 0), templateAffine)

    return morphed_mask_name



def rish_stat(siteName, imgs, templatePath, templateAffine, N_shm):

    for i in range(0, N_shm+1, 2):
        imgData= []
        for imgPath in imgs:
            prefix = os.path.basename(imgPath).split('.')[0]
            directory = os.path.dirname(imgPath)
            imgData.append(load_nifti(os.path.join(directory, 'harm', f'{prefix}_WarpedL{i}.nii.gz'))[0])

        save_nifti(os.path.join(templatePath, f'Mean_{siteName}_L{i}.nii.gz'),
                            np.mean(imgData, axis= 0), templateAffine)

        save_nifti(os.path.join(templatePath, f'Std_{siteName}_L{i}.nii.gz'),
                            np.std(imgData, axis= 0), templateAffine)


def template_masking(refMaskPath, targetMaskPath, templatePath, siteName, diffusionMeasures):

    refMask, affine= load_nifti(refMaskPath)
    targetMask, _= load_nifti(targetMaskPath)
    templateMask= applymask(refMask, targetMask)

    save_nifti(os.path.join(templatePath, 'templateMask.nii.gz'),
               templateMask, affine)

    for dm in diffusionMeasures:
        fileName= os.path.join(templatePath, f'Mean_{siteName}_{dm}.nii.gz')
        imgData, affine= load_nifti(fileName)
        save_nifti(fileName, applymask(imgData, templateMask), affine)

        fileName= os.path.join(templatePath, f'Std_{siteName}_{dm}.nii.gz')
        imgData, affine= load_nifti(fileName)
        save_nifti(fileName, applymask(imgData, templateMask), affine)


    return templateMask

def clip(data, L, H):
    data[np.isnan(data)]= 0.
    data= np.clip(data, L, H)

    return data

def smooth(data):

    # Equivalence b/w a Gaussian and a Box filter
    # https://stackoverflow.com/questions/35340197/box-filter-size-in-relation-to-gaussian-filter-sigma
    # Matlab does [3x3x3] box smoothing by default
    sigma= np.repeat(np.sqrt(2/3),3)

    return gaussian_filter(data, sigma)


def stat_calc(ref, target, mask):

    delta= applymask((ref- target), mask)
    per_diff= 100*delta/(ref+eps)
    per_diff= clip(per_diff, 100., -100.)
    per_diff_smooth= smooth(per_diff)
    scale= ref/(target+eps)

    return (delta, per_diff, per_diff_smooth, scale)


def difference_calc(refSite, targetSite, refImgs, targetImgs,
                    templatePath, templateAffine, subDir, mask, measures, travelHeads):

    '''
    if traveling heads:
        for each subject:
           delta= ref- target
        mean(delta)
        calc statistics
        save results
    else:
        delta= load(mean(ref))- load(mean(target))
        calc statistics
        save results
    '''

    for dm in measures:
        delta=[]
        per_diff=[]
        per_diff_smooth= []
        scale= []
        if travelHeads:
            for refImg, targetImg in zip(refImgs, targetImgs):
                prefix = os.path.basename(refImg).split('.')[0]
                directory = os.path.dirname(refImg)
                ref= load_nifti(os.path.join(directory, subDir, f'{prefix}_Warped{dm}.nii.gz'))[0]

                prefix = os.path.basename(targetImg).split('.')[0]
                directory = os.path.dirname(targetImg)
                target= load_nifti(os.path.join(directory, subDir, f'{prefix}_Warped{dm}.nii.gz'))[0]

                temp= stat_calc(ref, target, mask)
                delta.append(temp[0])
                per_diff.append(temp[1])
                per_diff_smooth.append(temp[2])
                scale.append(temp[3])

        else:
            fileName= os.path.join(templatePath, f'Mean_{refSite}_{dm}.nii.gz')
            ref= load_nifti(fileName)[0]

            fileName= os.path.join(templatePath, f'Mean_{targetSite}_{dm}.nii.gz')
            target= load_nifti(fileName)[0]

            temp= stat_calc(ref, target, mask)
            delta.append(temp[0])
            per_diff.append(temp[1])
            per_diff_smooth.append(temp[2])
            scale.append(temp[3])

        save_nifti(os.path.join(templatePath, f'Delta_{dm}.nii.gz'),
                   np.mean(delta, axis= 0), templateAffine)

        save_nifti(os.path.join(templatePath, f'PercentageDiff_{dm}.nii.gz'),
                   np.mean(per_diff, axis= 0), templateAffine)

        save_nifti(os.path.join(templatePath, f'PercentageDiff_{dm}smooth.nii.gz'),
                   np.mean(per_diff_smooth, axis= 0), templateAffine)

        save_nifti(os.path.join(templatePath, f'Scale_{dm}.nii.gz'),
                   np.sqrt(np.mean(scale, axis= 0)), templateAffine)



def main():
    import sys
    antsMult(sys.argv[1], sys.argv[2])

if __name__ == '__main__':
    main()

'''  
/home/tb571/Downloads/Harmonization-Python/lib/antsMult.py \
/home/tb571/Downloads/Harmonization-Python/connectom_prisma_demoData/template/CaseList_MultivariateTemplate.txt \
/home/tb571/Downloads/Harmonization-Python/lib/abc
'''