from subprocess import check_call
from buildTemplate import applyXform, modifiedFile
import os
from skimage.measure import label, regionprops
import numpy as np
from rish import rish

import warnings
with warnings.catch_warnings():
    warnings.filterwarnings("ignore", category=FutureWarning)

    from dipy.io.image import load_nifti, save_nifti
    from dipy.segment.mask import applymask

from scipy.ndimage import binary_erosion, binary_dilation, generate_binary_structure, median_filter

def antsReg(img, mask, mov, outPrefix):

    check_call((' ').join(['antsRegistrationSyNQuick.sh',
                           '-d', '3',
                           '-f', img,
                           '-x', mask,
                           '-m', mov,
                           '-o', outPrefix]), shell= True)

def antsApply(templatePath, directory, prefix, N_shm):

    for i in range(0, N_shm, 2):
        moving= os.path.join(templatePath, f'Scale_L{i}.nii.gz')
        # pass os.path.join(directory,harm)
        fixed= os.path.join(directory, 'harm', f'Scale_L{i}.nii.gz')
        output= os.path.join(directory, 'harm', f'Scale_L{i}_{prefix}.nii.gz')
        warp= os.path.join(directory, 'harm', f'ToSubjectSpace_{prefix}1Warp.nii.gz')
        trans= os.path.join(directory, 'harm', f'ToSubjectSpace_{prefix}0GenericAffine.mat')

        applyXform(moving, fixed, warp, trans, output)


def ring_masking(directory, prefix, maskPath, N_shm, shm_coeff, fit_matrix, b0):

    mapped_cs= []
    shs_same_level= [[0, 1], [1, 6], [6, 15], [15, 28], [28, 45]]

    for i in range(0, N_shm, 2):

        # load data and mask
        fileName= os.path.join(directory, f'Scale_L[{i}_{prefix}.nii.gz')
        img, affine = load_nifti(fileName)
        mask, _ = load_nifti(maskPath)

        if i==0: # compute the maskRing from 0th shm
            mask_scale= label(img>0.00001)
            maxArea= 0
            for region in regionprops(mask_scale):
                if region.area > maxArea:
                    maxLabel= region.label

            mask_scale= (mask_scale==maxLabel)
            mask*= mask_scale

            n_zero= 20
            se= generate_binary_structure(n_zero)
            maskTmp= np.pad(mask, (n_zero, n_zero, n_zero), 'constant', constant_values= 0.)

            dilM = binary_dilation(maskTmp, se)
            eroM = binary_erosion(maskTmp, se)
            maskRing = dilM - eroM

        img = applymask(img, mask)

        scaleTmp = np.pad(img, (n_zero, n_zero, n_zero), 'constant', constant_values= 0.)
        imgRing = applymask(scaleTmp, maskRing)

        percentile_mask= imgRing>=np.percentile(scaleTmp[maskRing>0], 95)
        roi= applymask(scaleTmp, percentile_mask)
        tmp= (median_filter(roi, (5,5,5)), percentile_mask)
        img= tmp[n_zero:-n_zero, n_zero:-n_zero, n_zero:-n_zero]

        save_nifti(fileName, data= img, affine= affine)

        ind= int(i/2)
        for level in shs_same_level[ind]:
            mapped_cs.append(img * shm_coeff[ :,:,:,level])


    S_hat= np.moveaxis(mapped_cs, 0, -1) @ fit_matrix
    S_hat[S_hat<0]= 0
    S_hat[S_hat>1]= 1

    # affine= templateAffine for all Scale_L{i}
    mappedFile= os.path.join(directory, 'mapped_cs.nii.gz')
    save_nifti(mappedFile, S_hat, affine= affine)

    # un-normalize harmonized data
    # load b0
    # multiply by b0
    # save result
    S_hat_dwi= applymask(S_hat, b0) # overriding applymask function with a nonbinary mask b0
    S_hat_final= np.concatenate(np.expand_dims(b0, axis=3), S_hat_dwi)
    mappedFile= os.path.join(directory, f'harmonized{prefix}.nii.gz')
    save_nifti(mappedFile, S_hat_final, affine= affine)

    return mappedFile

