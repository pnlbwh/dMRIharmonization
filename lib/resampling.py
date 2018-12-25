#!/usr/bin/env python

import warnings
with warnings.catch_warnings():
    warnings.filterwarnings("ignore", category=FutureWarning)

    from dipy.io.image import load_nifti, save_nifti
    from dipy.segment.mask import applymask
    import nibabel as nib

import numpy as np
from scipy.interpolate import splprep, splev, RegularGridInterpolator
from scipy.ndimage import binary_opening, generate_binary_structure
eps= 2.2204e-16

def RemoveGibbsRinging(b0):
    pass
    b0_n= None
    return b0_n

def save_high_res(fileName, sp_high, lowResImgHdr, highResImg):
    imgHdrOut = lowResImgHdr.copy()
    imgHdrOut['pixdim'][1:4] = imgHdrOut['pixdim'][1:4] * sp_high
    imgHdrOut['dim'][1:4] = highResImg.shape[1:4]
    save_nifti(fileName, highResImg, affine=imgHdrOut.get_qform(), header=imgHdrOut)

def resample(lowResImgPath, lowResMaskPath, bvals, sp_high, sOrder):

    # resample the dwi
    lowRes = nib.load(lowResImgPath)
    lowResImg= lowRes.get_data()
    lowResImgHdr= lowRes.hdr

    sp_low= lowResImgHdr['pixdim'][1:4]
    sx, sy, sz, n = lowResImg.shape
    step= sp_high/sp_low

    x = np.arange(0, sx, step[0])
    y = np.arange(0, sy, step[1])
    z = np.arange(0, sz, step[2])

    X, Y, Z= np.meshgrid(x,y,z)
    highResImg= np.zeros((X.shape[0], X.shape[1], X.shape[2], lowResImg.shape[3]), dtype='float')
    for i in range(lowResImg.shape[3]):
        img= lowResImg[:,:,:,i]
        tck = splprep([img[:, :, 0].flatten(), img[:, :, 1].flatten(), img[:, :, 2].flatten()], k=sOrder)[0]
        highResImg[:,:,:,i]= splev([X, Y, Z], tck).reshape(img.shape)

    # resample the mask
    lowRes = nib.load(lowResImgPath)
    lowResMask= lowRes.get_data()
    lowResMaskHdr= lowRes.hdr
    interp3= RegularGridInterpolator((np.arange(sx), np.arange(sy), np.arange(sz)), lowResMask)
    highResMask= interp3(([X, Y, Z]))
    highResMask= binary_opening(highResMask >= 0.5, structure=generate_binary_structure(3, 1)) * 1

    # resample the b0

    # TODO the following
    # extract b0
    where_b0= np.where(bvals == 0)[0]
    b0= lowResImg[...,where_b0].mean(-1)

    # process b0
    b0= None
    tck = splprep([b0[:, :, 0].flatten(), b0[:, :, 1].flatten(), b0[:, :, 2].flatten()], k=sOrder)[0]
    b0_n= splev([X, Y, Z], tck).reshape(b0.shape)
    b0_n= b0_n*highResMask

    b0_n[b0_n<=0]= 1.
    b0_n[np.isnan(b0_n)]= 1.

    # pass the following step due to unavailabily of Python code
    # b0_gibs = RemoveGibbsRinging(b0_n)
    b0_gibs= b0_n

    ind= np.where(lowResImg>0)[0]
    b0_gibs[b0_gibs > b0.max()] = b0[ind].max()
    b0_gibs[b0_gibs < b0.min()] = b0[ind].min()

    b0_gibs[b0_gibs <= 0] = 1.
    b0_gibs[np.isnan(b0_gibs)] = 1.

    # insert b0 back
    highResImg[:,:,:,where_b0]= b0_gibs

    highResImg[highResImg > (lowResImg+eps).max()] = lowResImg[ind,:].max()
    highResImg[highResImg < (lowResImg-eps).min()] = lowResImg[ind,:].min()

    highResImg= applymask(highResImg, highResMask)

    highResImg[highResImg < eps] = 0.
    highResImg[highResImg > 1] = 1.

    highResImgPath= lowResImgPath.split('.')[0]+'_resampled'+'.nii.gz'
    highResMaskPath= lowResMaskPath.split('.')[0]+'_resampled'+'.nii.gz'

    save_high_res(highResImgPath, sp_high, lowResImgHdr, highResImg)
    save_high_res(highResMaskPath, sp_high, lowResMaskHdr, highResMask)
    save_high_res(highResImgPath.split('.')[0]+'_b0'+'.nii.gz', sp_high, lowResMaskHdr, b0_gibs)

    return (highResImgPath, highResMaskPath)


if __name__=='__main__':
    lowResImgPath= '/home/tb571/Downloads/Harmonization-Python/connectom_prisma_demoData/A/connectom/dwi_A_connectom_st_b1200.nii.gz'
    lowResMaskPath= '/home/tb571/Downloads/Harmonization-Python/connectom_prisma_demoData/A/connectom/mask.nii.gz'

    from dipy.io import read_bvals_bvecs
    bvals, _= read_bvals_bvecs(
        '/home/tb571/Downloads/Harmonization-Python/connectom_prisma_demoData/A/connectom/dwi_A_connectom_st_b1200.bval',
        '/home/tb571/Downloads/Harmonization-Python/connectom_prisma_demoData/A/connectom/dwi_A_connectom_st_b1200.bvec')

    resample(lowResImgPath, lowResMaskPath, bvals, np.array([1.7188,1.7187,3]), 5)

