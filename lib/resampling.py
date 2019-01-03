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
from subprocess import check_call
eps= 2.2204e-16


def Bspline(lowFile, highFile, sp_high):
    check_call((' ').join(['ResampleImage', '3',
                           lowFile, highFile,
                           f'{sp_high[0]}x{sp_high[1]}x{sp_high[2]}',
                           '1', '4', '5']), shell= True)

    return nib.load(highFile).get_data()

def Linear(lowFile, highFile, sp_high):
    check_call((' ').join(['ResampleImage', '3',
                           lowFile, highFile,
                           f'{sp_high[0]}x{sp_high[1]}x{sp_high[2]}',
                           '1']), shell= True)

    return nib.load(highFile).get_data()

def RemoveGibbsRinging(lowResB0, higResB0, sp_high):

    # resample the b0
    Bspline(lowResB0, higResB0, sp_high)

    # unrign the b0
    check_call(['unring.a64', higResB0, higResB0])

    # load b0 back
    return nib.load(higResB0).get_data()


def save_high_res(fileName, sp_high, lowResImgHdr, highResImg):
    imgHdrOut = lowResImgHdr.copy()
    imgHdrOut['pixdim'][1:4] = sp_high
    imgHdrOut['dim'][1:4] = highResImg.shape[:3]
    save_nifti(fileName, highResImg, affine=imgHdrOut.get_qform(), hdr=imgHdrOut)

def resample(lowResImgPath, lowResMaskPath, bvals, sp_high, sOrder):

    # resample the dwi
    lowRes = nib.load(lowResImgPath)
    lowResImg= lowRes.get_data()
    lowResImgHdr= lowRes.header

    sp_low= lowResImgHdr['pixdim'][1:4]
    sx, sy, sz, n = lowResImg.shape
    step= sp_high/sp_low

    x = np.arange(0, sx)
    y = np.arange(0, sy)
    z = np.arange(0, sz)
    X_old, Y_old, Z_old = np.meshgrid(x, y, z)

    x = np.arange(0, sx-1, step[0])
    y = np.arange(0, sy-1, step[1])
    z = np.arange(0, sz-1, step[2])
    X, Y, Z= np.meshgrid(x, y, z)

    highResImg= np.zeros((X.shape[0], X.shape[1], X.shape[2], lowResImg.shape[3]), dtype='float')
    for i in range(lowResImg.shape[3]):
        print(i)
        img= lowResImg[:,:,:,i]
        # tck = splprep([X_old.flatten(), Y_old.flatten(), Z_old.flatten(), img.flatten()], k=sOrder)[0]
        # highResImg[:,:,:,i]= splev([X.flatten(), Y.flatten(), Z.flatten()], tck).reshape(X.shape)
        # tmpFile= lowResImgPath.split('.')[0]+ 'tmp.nii.gz'
        tmpFile= 'tmp.nii.gz'
        save_high_res(tmpFile, sp_low, lowResImgHdr, img)
        highResImg[:, :, :, i]= Bspline(tmpFile, tmpFile, X.shape)

    # resample the mask
    lowRes = nib.load(lowResMaskPath)
    lowResMask= lowRes.get_data()
    lowResMaskHdr= lowRes.header
    # interp3= RegularGridInterpolator((np.arange(sx), np.arange(sy), np.arange(sz)), lowResMask)
    # temp= []
    # for xi in x:
    #     for yi in y:
    #         for zi in z:
    #             print(xi, yi, zi)
    #             temp.append(interp3([xi, yi, zi]))
    # highResMask = np.reshape(temp, X.shape)

    highResMaskPath = lowResMaskPath.split('.')[0] + '_resampled' + '.nii.gz'
    highResMask= Linear(lowResMaskPath, highResMaskPath, X.shape)

    highResMask= binary_opening(highResMask >= 0.5, structure=generate_binary_structure(3, 1)) * 1


    save_high_res(highResMaskPath, sp_high, lowResMaskHdr, highResMask)


    # resample the b0


    # extract b0
    where_b0= np.where(bvals == 0)[0]
    b0= lowResImg[...,where_b0].mean(-1)

    # process b0
    # tck = splprep([b0[:, :, 0].flatten(), b0[:, :, 1].flatten(), b0[:, :, 2].flatten()], k=sOrder)[0]
    # b0_n= splev([X, Y, Z], tck).reshape(b0.shape)
    # b0_n= b0_n*highResMask

    # b0_n[b0_n<=0]= 1.
    # b0_n[np.isnan(b0_n)]= 1.

    # pass the following step due to unavailabily of Python code
    lowResB0 = lowResImgPath.split('.')[0] + '_bse' + '.nii.gz'
    highResB0= lowResImgPath.split('.')[0] + '_bse_resampled' + '.nii.gz'
    b0_gibs = RemoveGibbsRinging(lowResB0, highResB0, X.shape)

    ind_x, ind_y, ind_z, _= np.where(lowResImg>0)
    b0_gibs[b0_gibs > b0.max()] = b0[ind_x, ind_y, ind_z].max()
    b0_gibs[b0_gibs < b0.min()] = b0[ind_x, ind_y, ind_z].min()

    # b0_gibs[b0_gibs <= 0] = 1.
    # b0_gibs[np.isnan(b0_gibs)] = 1.
    save_high_res(highResB0, sp_high, lowResMaskHdr, b0_gibs)

    # insert b0 back
    for i in where_b0:
        highResImg[:,:,:,i]= b0_gibs

    highResImg[highResImg > (lowResImg+eps).max()] = lowResImg[ind_x, ind_y, ind_z,: ].max()
    highResImg[highResImg < (lowResImg-eps).min()] = lowResImg[ind_x, ind_y, ind_z,: ].min()

    highResImg= applymask(highResImg, highResMask)

    # highResImg[highResImg < eps] = 0.
    # highResImg[highResImg > 1] = 1.

    highResImgPath= lowResImgPath.split('.')[0]+'_resampled'+'.nii.gz'
    save_high_res(highResImgPath, sp_high, lowResImgHdr, highResImg)

    return (highResImgPath, highResMaskPath)


if __name__=='__main__':
    from dipy.io import read_bvals_bvecs

    # lowResImgPath= '/home/tb571/Downloads/Harmonization-Python/connectom_prisma_demoData/A/connectom/dwi_A_connectom_st_b1200.nii.gz'
    # lowResMaskPath= '/home/tb571/Downloads/Harmonization-Python/connectom_prisma_demoData/A/connectom/mask.nii.gz'
    # bvals, _= read_bvals_bvecs(
    #     '/home/tb571/Downloads/Harmonization-Python/connectom_prisma_demoData/A/connectom/dwi_A_connectom_st_b1200.bval',
    #     '/home/tb571/Downloads/Harmonization-Python/connectom_prisma_demoData/A/connectom/dwi_A_connectom_st_b1200.bvec')

    lowResImgPath= '/home/tb571/Downloads/Harmonization-Python/test_data/test_a/connectom/connectom_a_dwi.nii.gz'
    lowResMaskPath= '/home/tb571/Downloads/Harmonization-Python/test_data/test_a/connectom/connectom_a_mask.nii.gz'
    bvals, _= read_bvals_bvecs(
        '/home/tb571/Downloads/Harmonization-Python/test_data/test_a/connectom/connectom_a_dwi.bval',
        '/home/tb571/Downloads/Harmonization-Python/test_data/test_a/connectom/connectom_a_dwi.bvec')


    resample(lowResImgPath, lowResMaskPath, bvals, np.array([1.7188,1.7187,3]), 7)

