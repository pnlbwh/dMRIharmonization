#!/usr/bin/env python

import numpy as np
from subprocess import check_call
from plumbum.cmd import antsApplyTransforms
from plumbum.cmd import FG
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

    # dir_bak = os.getcwd()

    prefix= os.path.basename(imgPath).split('.')[0]
    warp = glob(os.path.join(templatePath, prefix + f'_{dm}*[!Inverse]Warp.nii.gz'))
    trans = glob(os.path.join(templatePath, prefix + f'_{dm}*GenericAffine.mat'))

    # warping the mask
    applyXform(maskPath,
               os.path.join(templatePath, 'template0.nii.gz'),
               warp, trans,
               str(maskPath).split('.')[0]+ 'Warped.nii.gz')


    # warping the rish features
    # os.chdir(rishPath)
    for i in range(0, N_shm, 2):
        applyXform(modifiedFile(imgPath, 'harm', f'_L{i}.nii.gz'),
                   os.path.join(templatePath, 'template0.nii.gz'),
                   warp, trans,
                   modifiedFile(imgPath, 'harm', f'_WarpedL{i}.nii.gz'))



    # warping the diffusion measures
    # os.chdir(dtiPath)
    for dm in diffusionMeasures:
        applyXform(modifiedFile(imgPath, 'dti', f'_{dm}.nii.gz'),
                   os.path.join(templatePath, 'template0.nii.gz'),
                   warp, trans,
                   modifiedFile(imgPath, 'dti', f'_Warped{dm}.nii.gz'))


    # os.chdir(dir_bak)


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


def modifiedFile(file, subDir, ext):
    prefix= os.path.basename(file).split('.')[0]
    return os.path.join(os.path.dirname(file), subDir, prefix+ext)


def dti_stat(siteName, imgs, masks, templatePath, templateAffine, diffusionMeasures):

    maskData = []
    for mask in masks:
        maskData.append(load_nifti(modifiedFile(mask, 'harm', 'Warped.nii.gz'))[0])

    morphed_mask= binary_opening(np.mean(maskData, axis=0)>0.5, structure= generate_binary_structure(3,1))
    morphed_mask_name= os.path.join(templatePath, f'{siteName}_Mask.nii.gz')
    save_nifti(morphed_mask_name, morphed_mask, templateAffine, None)

    imgData= []
    for dm in diffusionMeasures:
        for img in imgs:
            imgData.append(load_nifti(modifiedFile(img, 'dti', f'_Warped{dm}.nii.gz'))[0])

        save_nifti(os.path.join(templatePath, f'Mean_{siteName}_{dm}.nii.gz'),
                                np.mean(imgData), templateAffine, None)

        save_nifti(os.path.join(templatePath, f'Std_{siteName}_{dm}.nii.gz'),
                                np.std(imgData), templateAffine, None)

    return morphed_mask_name



# extns= diffusionMeasures
# extns= [f'L{i}' for i in range(0, N_shm, 2)]

# def mean_std_calc(imgs, subDir, extns, templatePath):
#     imgData= []
#     for dm in diffusionMeasures:
#         for img in imgs:
#             imgData.append(load_nifti(modifiedFile(img, subDir, f'_Warped{ext}.nii.gz'))[0])
#
#         save_nifti(os.path.join(templatePath, f'Mean_{siteName}_Warped{dm}.nii.gz'),
#                                 np.mean(imgData), templateAffine, None)
#
#         save_nifti(os.path.join(templatePath, f'Std_{siteName}_Warped{dm}.nii.gz'),
#                                 np.std(imgData), templateAffine, None)


def rish_stat(siteName, imgs, templatePath, templateAffine, N_shm):

    for i in range(0, N_shm, 2):
        imgData= []
        for img in imgs:
            imgData.append(load_nifti(modifiedFile(img, 'harm', f'_WarpedL{i}.nii.gz'))[0])

        save_nifti(os.path.join(templatePath, f'Mean_{siteName}_L{i}.nii.gz'),
                            np.mean(imgData), templateAffine, None)

        save_nifti(os.path.join(templatePath, f'Std_{siteName}_L{i}.nii.gz'),
                            np.std(imgData), templateAffine, None)


def template_masking(refMaskPath, targetMaskPath, templatePath, siteName, diffusionMeasures):

    refMask, affine= load_nifti(refMaskPath)
    targetMask, _= load_nifti(targetMaskPath)
    templateMask= applymask(refMask, targetMask)

    save_nifti(os.path.join(templatePath, 'templateMask.nii.gz'),
               templateMask, affine)

    for dm in diffusionMeasures:
        fileName= os.path.join(templatePath, f'Mean_{siteName}_Warped{dm}.nii.gz')
        imgData, affine= load_nifti(fileName)
        save_nifti(fileName, applymask(imgData, templateMask), affine, None)

        fileName= os.path.join(templatePath, f'Std_{siteName}_Warped{dm}.nii.gz')
        imgData, affine= load_nifti(fileName)
        save_nifti(fileName, applymask(imgData, templateMask), affine, None)


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


# call once for dti:
# difference_calc('reference', 'target', refImgs, targetImgs, templatePath, templateAffine,
#                 'dti', templateMask, diffusionMeasures, travelHeads)

# call again for rish:
# difference_calc('reference', 'target', refImgs, targetImgs, templatePath, templateAffine,
#                 'harm', templateMask, [f'L{i}' for i in range(0, N_shm, 2)], travelHeads)

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
        delta=[], per_diff=[], per_diff_smooth= [], scale= []
        if travelHeads:
            for refImg, targetImg in (refImgs, targetImgs):
                ref= load_nifti(modifiedFile(refImg, subDir, f'_Warped{dm}.nii.gz'))[0]
                target= load_nifti(modifiedFile(targetImg, subDir, f'_Warped{dm}.nii.gz'))[0]

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
                   delta, templateAffine)

        save_nifti(os.path.join(templatePath, f'PercentageDiff_{dm}.nii.gz'),
                   per_diff, templateAffine)

        save_nifti(os.path.join(templatePath, f'PercentageDiff_{dm}smooth.nii.gz'),
                   per_diff_smooth, templateAffine)

        save_nifti(os.path.join(templatePath, f'Scale_{dm}.nii.gz'),
                   scale, templateAffine)



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