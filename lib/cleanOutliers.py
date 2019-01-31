from subprocess import check_call
from buildTemplate import applyXform
import os
from skimage.measure import label, regionprops
from rish import rish
from local_med_filter import local_med_filter

import warnings
with warnings.catch_warnings():
    warnings.filterwarnings("ignore", category=FutureWarning)

    from dipy.io.image import load_nifti, save_nifti
    from dipy.segment.mask import applymask

import numpy as np

from scipy.ndimage import binary_erosion, binary_dilation, \
    generate_binary_structure, median_filter, iterate_structure
eps= 2.2204e-16

def antsReg(img, mask, mov, outPrefix):

    if mask:
        check_call((' ').join(['antsRegistrationSyNQuick.sh',
                               '-d', '3',
                               '-f', img,
                               '-x', mask,
                               '-m', mov,
                               '-o', outPrefix]), shell= True)
    else:
        check_call((' ').join(['antsRegistrationSyNQuick.sh',
                               '-d', '3',
                               '-f', img,
                               '-m', mov,
                               '-o', outPrefix]), shell= True)

def antsApply(templatePath, directory, prefix, N_shm):

    for i in range(0, N_shm+1, 2):
        moving= os.path.join(templatePath, f'Scale_L{i}.nii.gz')
        fixed= os.path.join(directory, f'{prefix}_L{i}.nii.gz')
        output= os.path.join(directory, f'Scale_L{i}_{prefix}.nii.gz')
        warp= os.path.join(directory, f'ToSubjectSpace_{prefix}1Warp.nii.gz')
        trans= os.path.join(directory, f'ToSubjectSpace_{prefix}0GenericAffine.mat')

        applyXform(moving, fixed, warp, trans, output)


def custom_spherical_structure(D):
    '''Given diameter D, returns a spherical structuring element'''

    sw=(D-1)/2
    ses2= int(np.ceil(D/2))
    [y,x,z]= np.meshgrid(np.arange(-sw,sw+1),np.arange(-sw,sw+1),np.arange(-sw,sw+1))
    m= np.sqrt(x**2 + y**2 + z**2)
    b= (m <= m[ses2,ses2,D-1])

    return b


def findLargestConnectMask(img, mask):

    mask_scale = label(img > 0.00001, connectivity=1)
    maxArea = 0
    for region in regionprops(mask_scale):
        if region.area > maxArea:
            maxLabel = region.label
            maxArea = region.area

    largeConnectMask = (mask_scale == maxLabel)
    mask *= largeConnectMask

    return mask


def ring_masking(directory, prefix, maskPath, N_shm, shm_coeff, b0, qb_model):

    B = qb_model.B
    bvals= qb_model.gtab.bvals

    mapped_cs= []
    shs_same_level= [[0, 1], [1, 6], [6, 15], [15, 28], [28, 45]]

    for i in range(0, N_shm+1, 2):

        # load data and mask
        fileName= os.path.join(directory, 'harm', f'Scale_L{i}_{prefix}.nii.gz')
        img, affine = load_nifti(fileName)
        mask, _ = load_nifti(maskPath)


        if i==0: # compute the skullRingMask from 0th shm

            mask= findLargestConnectMask(img, mask)

            n_zero= 5
            se= custom_spherical_structure(n_zero)
            paddedMask= np.pad(mask, n_zero, 'constant', constant_values= 0.)

            dilM = binary_dilation(paddedMask, se)*1
            eroM = binary_erosion(paddedMask, se)*1
            skullRingMask = dilM - eroM

        paddedImg = np.pad(img, n_zero, 'constant', constant_values=0.)
        skullRing= paddedImg*skullRingMask
        thresh= np.percentile(skullRing[skullRing>0], 95)
        outLier= (skullRing>thresh)*1
        tmp= local_med_filter(paddedImg, outLier)
        denoisedImg = tmp[n_zero:-n_zero, n_zero:-n_zero, n_zero:-n_zero]

        # restImg= paddedImg*(1-skullRingMask)
        # denoisedSkullRing= median_filter(skullRing, (5,5,5))
        # # multiply denoisedSkullRing by skullRing again to avoid out of ROI inclusion
        # tmp= restImg+denoisedSkullRing*skullRing
        # denoisedImg= tmp[n_zero:-n_zero, n_zero:-n_zero, n_zero:-n_zero]


        # for i=0, create new mask over denoisedImg, save it as harmonzied_{prefix}_mask
        if i==0:
            mask_final= findLargestConnectMask(denoisedImg, mask)
            mask_final= binary_dilation(mask_final, generate_binary_structure(3,1))*1
            harmMask = os.path.join(directory, f'harmonized_{prefix}_mask.nii.gz')
            save_nifti(harmMask, mask_final, affine=affine)


        ind= int(i/2)
        for level in range(shs_same_level[ind][0], shs_same_level[ind][1]):
            mapped_cs.append(denoisedImg * shm_coeff[ :,:,:,level])


    # FIXME: incorrect product
    # S_hat= np.moveaxis(mapped_cs, 0, -1) @ B.T
    S_hat= np.dot(np.moveaxis(mapped_cs, 0, -1), B.T)
    S_hat[S_hat<eps]= 0
    S_hat[S_hat>1]= 1

    # affine= templateAffine for all Scale_L{i}
    mappedFile= os.path.join(directory, f'{prefix}_mapped_cs.nii.gz')
    save_nifti(mappedFile, S_hat, affine= affine)

    # un-normalize harmonized data
    S_hat_dwi= applymask(S_hat, b0) # overriding applymask function with a nonbinary mask b0

    # place b0s in proper indices
    S_hat_final= stack_b0(bvals, S_hat_dwi, b0)
    S_hat_final= applymask(S_hat_final, mask_final)
    # S_hat_final= np.concatenate((np.expand_dims(b0, axis=3), S_hat_dwi), axis= 3)


    # save harmonized data
    harmImg= os.path.join(directory, f'harmonized_{prefix}.nii.gz')
    save_nifti(harmImg, S_hat_final, affine= affine)


    return (harmImg, harmMask)


def stack_b0(bvals, dwi, b0):

    ind= np.where(bvals==0)[0]

    S_hat_final= []
    j= 0
    for i in range(len(bvals)):
        if i in ind:
            S_hat_final.append(b0)
        else:
            S_hat_final.append(dwi[:,:,:,j])
            j+=1

    return np.moveaxis(np.array(S_hat_final), 0, -1)

if __name__ =='__main__':
    imgPath= '/home/tb571/Downloads/Harmonization-Python/connectom_prisma_demoData/A/prisma/dwi_A_prisma_st_b1200.nii.gz'
    maskPath= '/home/tb571/Downloads/Harmonization-Python/connectom_prisma_demoData/A/prisma/mask.nii.gz'

    directory= os.path.dirname(imgPath)
    inPrefix= imgPath.split('.')[0]
    prefix= os.path.split(inPrefix)[-1]
    outPrefix = os.path.join(directory, 'harm', prefix)
    N_shm= 6

    b0, shm_coeff, qb_model= rish(imgPath, maskPath, inPrefix, outPrefix, N_shm)
    ring_masking(directory, prefix, maskPath, N_shm, shm_coeff, b0, qb_model)