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

from skimage.transform import resize
from scipy.ndimage import binary_opening, generate_binary_structure
from normalize import normalize_data, find_b0
from util import *

def save_high_res(fileName, sp_high, lowResImgHdr, highResImg):

    imgHdrOut = lowResImgHdr.copy()
    sp_low= imgHdrOut['pixdim'][1:4]
    imgHdrOut['pixdim'][1:4] = sp_high
    imgHdrOut['dim'][1:4] = highResImg.shape[:3]
    scale= np.diag((sp_high/sp_low).tolist()+[1.])
    imgHdrOut.set_sform(imgHdrOut.get_sform() @ scale)
    imgHdrOut.set_qform(imgHdrOut.get_qform() @ scale)
    save_nifti(fileName, highResImg, affine= imgHdrOut.get_best_affine(), hdr=imgHdrOut)


def resampling(lowResImgPath, lowResMaskPath, lowResImg, lowResImgHdr, lowResMask, lowResMaskHdr, sp_high, bvals):

    # order for b spline interpolation
    sOrder= 5

    # resample the dwi ----------------------------------------------------------------
    where_b0 = np.where(bvals == 0)[0]
    lowResImgPrime, b0= normalize_data(lowResImg, mask= lowResMask, where_b0= where_b0)
    lowResImg= applymask(lowResImgPrime, b0)

    sp_low= lowResImgHdr['pixdim'][1:4]
    step = sp_low/sp_high
    sx, sy, sz = [int(round(x)) for x in lowResImg.shape[:3]*step]

    highResImg= np.zeros((sx, sy, sz, lowResImg.shape[3]), dtype='float')
    for i in np.where(bvals != 0)[0]:
        print('Resampling gradient ', i)
        highResImg[:,:,:,i]= resize(lowResImg[:,:,:,i], (sx, sy, sz), order= sOrder, mode= 'constant')

    # resample the mask ---------------------------------------------------------------
    highResMaskPath = lowResMaskPath.split('.')[0] + '_resampled.nii.gz'
    highResMask= resize(lowResMask.astype('float'), (sx, sy, sz), order= 1, mode= 'constant') # order 1 for linear interpolation
    highResMask= binary_opening(highResMask >= 0.5, structure=generate_binary_structure(3, 1)) * 1
    save_high_res(highResMaskPath, sp_high, lowResMaskHdr, highResMask.astype('uint8'))


    # resample the b0 ----------------------------------------------------------------
    highResB0PathTmp= lowResImgPath.split('.')[0] + '_resampled_bse_tmp.nii.gz'
    b0HighRes= resize(b0, (sx, sy, sz), order= sOrder, mode= 'constant')
    np.nan_to_num(b0HighRes).clip(min= 0., out= b0HighRes) # using min= 1. is unnecessary
    save_high_res(highResB0PathTmp, sp_high, lowResMaskHdr, b0HighRes)

    # unring the b0
    highResB0Path = lowResImgPath.split('.')[0] + '_resampled_bse.nii.gz'
    check_call(['unring.a64', highResB0PathTmp, highResB0Path])
    check_call(['rm', highResB0PathTmp])
    b0_gibs = load(highResB0Path).get_data()
    np.nan_to_num(b0_gibs).clip(min= 0., out= b0_gibs) # using min= 1. is unnecessary

    # defining lh_max and lh_min separately to deal with memory error
    lh_max= b0.max()
    lh_min= b0.min()
    b0_gibs[b0_gibs > lh_max] = lh_max
    b0_gibs[b0_gibs < lh_min] = lh_min
    save_high_res(highResB0Path, sp_high, lowResMaskHdr, b0_gibs)


    # insert b0 back ------------------------------------------------------------------
    for i in where_b0:
        highResImg[:,:,:,i]= b0_gibs

    lh_max= lowResImg.max()
    lh_min= lowResImg.min()
    highResImg[highResImg > lh_max] = lh_max
    highResImg[highResImg < lh_min] = lh_min

    highResImg= applymask(highResImg, highResMask)
    highResImgPath= lowResImgPath.split('.')[0]+'_resampled.nii.gz'
    save_high_res(highResImgPath, sp_high, lowResImgHdr, highResImg)

    return (highResImgPath, highResMaskPath)


if __name__=='__main__':
    pass

