#!/usr/bin/env python

import warnings
with warnings.catch_warnings():
    warnings.filterwarnings("ignore", category=FutureWarning)

    from dipy.io.image import load_nifti, save_nifti
    # from dipy.reconst.shm import normalize_data
    from dipy.segment.mask import applymask
    import nibabel as nib

import numpy as np
from skimage.transform import resize
from scipy.ndimage import binary_opening, generate_binary_structure
from subprocess import check_call
from normalize import normalize_data, find_b0
eps= 2.2204e-16


def save_high_res(fileName, sp_high, lowResImgHdr, highResImg):

    imgHdrOut = lowResImgHdr.copy()
    sp_low= imgHdrOut['pixdim'][1:4]
    imgHdrOut['pixdim'][1:4] = sp_high
    imgHdrOut['dim'][1:4] = highResImg.shape[:3]
    scale= np.diag((sp_high/sp_low).tolist()+[1.])
    imgHdrOut.set_sform(imgHdrOut.get_sform() @ scale)
    imgHdrOut.set_qform(imgHdrOut.get_qform() @ scale)
    save_nifti(fileName, highResImg, affine= imgHdrOut.get_qform(), hdr=imgHdrOut)

def resampling(lowResImgPath, lowResMaskPath, lowResImg, lowResImgHdr, lowResMask, lowResMaskHdr, sp_high, bvals):

    # order for b spline interpolation
    sOrder= 5

    # resample the dwi ----------------------------------------------------------------
    lowResImg= applymask(lowResImg, lowResMask)
    where_b0= np.where(bvals == 0)[0]
    b0= find_b0(lowResImg, where_b0= where_b0)

    sp_low= lowResImgHdr['pixdim'][1:4]
    step = sp_low/sp_high
    sx, sy, sz = [int(round(x)) for x in lowResImg.shape[:3]*step]

    highResImg= np.zeros((sx, sy, sz, lowResImg.shape[3]), dtype='float')
    for i in np.where(bvals != 0)[0]:
        print('Resampling gradient ', i)
        highResImg[:,:,:,i]= resize(np.double(lowResImg[:,:,:,i]), (sx, sy, sz), order= sOrder, mode= 'constant')

    # resample the mask ---------------------------------------------------------------
    highResMaskPath = lowResMaskPath.split('.')[0] + '_resampled' + '.nii.gz'
    highResMask= resize(np.double(lowResMask), (sx, sy, sz), order= 1, mode= 'constant') # order 1 for linear interpolation
    highResMask= binary_opening(highResMask >= 0.5, structure=generate_binary_structure(3, 1)) * 1
    save_high_res(highResMaskPath, sp_high, lowResMaskHdr, highResMask.astype(int))


    # resample the b0 ----------------------------------------------------------------
    highResB0Path= lowResImgPath.split('.')[0] + '_resampled_bse' + '.nii.gz'
    b0HighRes= resize(np.double(b0), (sx, sy, sz), order= sOrder, mode= 'constant')
    np.nan_to_num(b0HighRes).clip(min= 1., out= b0HighRes)
    save_high_res(highResB0Path, sp_high, lowResMaskHdr, b0HighRes)

    # unring the b0
    check_call(['unring.a64', highResB0Path, highResB0Path])
    b0_gibs = nib.load(highResB0Path).get_data()

    # defining lh_max and lh_min separately to deal with memory error
    lh_max= b0.max()
    lh_min= b0.min()
    b0_gibs[b0_gibs > lh_max] = lh_max
    b0_gibs[b0_gibs < lh_min] = lh_min
    np.nan_to_num(b0_gibs).clip(min= 1., out= b0_gibs)
    save_high_res(highResB0Path, sp_high, lowResMaskHdr, b0_gibs)


    # insert b0 back ------------------------------------------------------------------
    for i in where_b0:
        highResImg[:,:,:,i]= b0_gibs

    lh_max= lowResImg.max()
    lh_min= lowResImg.min()
    highResImg[highResImg > lh_max] = lh_max
    highResImg[highResImg < lh_min] = lh_min

    highResImg= applymask(highResImg, highResMask)
    highResImgPath= lowResImgPath.split('.')[0]+'_resampled'+'.nii.gz'
    highResImg, _= normalize_data(highResImg, b0= b0_gibs)
    highResImg= applymask(highResImg, b0_gibs)
    save_high_res(highResImgPath, sp_high, lowResImgHdr, highResImg)

    return (highResImgPath, highResMaskPath)


if __name__=='__main__':
    from dipy.io import read_bvals_bvecs

    # lowResImgPath= '/home/tb571/Downloads/Harmonization-Python/connectom_prisma_demoData/A/connectom/dwi_A_connectom_st_b1200.nii.gz'
    # lowResMaskPath= '/home/tb571/Downloads/Harmonization-Python/connectom_prisma_demoData/A/connectom/mask.nii.gz'
    # bvals, _= read_bvals_bvecs(
    #     '/home/tb571/Downloads/Harmonization-Python/connectom_prisma_demoData/A/connectom/dwi_A_connectom_st_b1200.bval',
    #     None)

    # lowResImgPath= '/home/tb571/Downloads/Harmonization-Python/test_data/test_a/connectom/connectom_a_dwi.nii.gz'
    # lowResMaskPath= '/home/tb571/Downloads/Harmonization-Python/test_data/test_a/connectom/connectom_a_mask.nii.gz'
    # bvals, _= read_bvals_bvecs(
    #     '/home/tb571/Downloads/Harmonization-Python/test_data/test_a/connectom/connectom_a_dwi.bval', None)

    lowResImgPath= '/home/tb571/Downloads/Harmonization-Python/BSNIP_Baltimore/BSNIP_Balt_trainingHC/GT_3507/GT_3507_dwi_xc_Ed.nii.gz'
    lowResMaskPath= '/home/tb571/Downloads/Harmonization-Python/BSNIP_Baltimore/BSNIP_Balt_trainingHC/GT_3507/GT_3507_dwi_xc_Ed_OTSUtensormask_cleaned.nii.gz'
    bvals, _= read_bvals_bvecs(
        '/home/tb571/Downloads/Harmonization-Python/BSNIP_Baltimore/BSNIP_Balt_trainingHC/GT_3507/GT_3507_dwi_xc_Ed.bval', None)

    lowRes = nib.load(lowResImgPath)
    lowResImg= lowRes.get_data()
    lowResImgHdr= lowRes.header

    lowRes = nib.load(lowResMaskPath)
    lowResMask= lowRes.get_data()
    lowResMaskHdr= lowRes.header

    resampling(lowResImgPath, lowResMaskPath, lowResImg, lowResImgHdr, lowResMask, lowResMaskHdr,
               np.array([1.5,1.5,1.5]), bvals)


