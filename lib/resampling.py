#!/usr/bin/env python

import warnings
with warnings.catch_warnings():
    warnings.filterwarnings("ignore", category=FutureWarning)

    from dipy.io.image import load_nifti, save_nifti
    from dipy.reconst.shm import normalize_data
    from dipy.segment.mask import applymask
    import nibabel as nib

import numpy as np
from skimage.transform import resize
from scipy.ndimage import binary_opening, generate_binary_structure
from subprocess import check_call
eps= 2.2204e-16


# unused
def Bspline(lowFile, highFile, sp_high):
    check_call((' ').join(['ResampleImage', '3',
                           lowFile, highFile,
                           f'{sp_high[0]}x{sp_high[1]}x{sp_high[2]}',
                           '1', '4', '5']), shell= True)

    return nib.load(highFile).get_data()

# unused
def Linear(lowFile, highFile, sp_high):
    check_call((' ').join(['ResampleImage', '3',
                           lowFile, highFile,
                           f'{sp_high[0]}x{sp_high[1]}x{sp_high[2]}',
                           '1']), shell= True)

    return nib.load(highFile).get_data()


def save_high_res(fileName, sp_high, lowResImgHdr, highResImg):
    imgHdrOut = lowResImgHdr.copy()
    imgHdrOut['pixdim'][1:4] = sp_high
    imgHdrOut['dim'][1:4] = highResImg.shape[:3]
    save_nifti(fileName, highResImg, affine=imgHdrOut.get_qform(), hdr=imgHdrOut)

def resampling(lowResImgPath, lowResMaskPath, lowResImg, lowResImgHdr, lowResMask, lowResMaskHdr, sp_high, bvals):

    # order for b spline interpolation
    sOrder= 5

    # resample the dwi ----------------------------------------------------------------
    lowResImg= applymask(lowResImg, lowResMask)
    where_b0= np.where(bvals == 0)[0]
    # b0_orig = lowResImg[..., where_b0].mean(-1)
    b0 = lowResImg[..., where_b0].mean(-1)


    sp_low= lowResImgHdr['pixdim'][1:4]
    step = sp_low/sp_high
    sx, sy, sz = [int(round(x)) for x in lowResImg.shape[:3]*step]

    highResImg= np.zeros((sx, sy, sz, lowResImg.shape[3]), dtype='float')
    for i in range(lowResImg.shape[3]):
        print('Resampling gradient ', i)
        highResImg[:,:,:,i]= resize(np.double(lowResImg[:,:,:,i]), (sx, sy, sz), order= sOrder)

    # resample the mask ---------------------------------------------------------------
    highResMaskPath = lowResMaskPath.split('.')[0] + '_resampled' + '.nii.gz'
    highResMask= resize(np.double(lowResMask), (sx, sy, sz), order= sOrder)
    highResMask= binary_opening(highResMask >= 0.5, structure=generate_binary_structure(3, 1)) * 1

    save_high_res(highResMaskPath, sp_high, lowResMaskHdr, highResMask.astype(int))


    # resample the b0 ----------------------------------------------------------------
    highResB0Path= lowResImgPath.split('.')[0] + '_resampled_bse' + '.nii.gz'
    b0HighRes= resize(np.double(b0), (sx, sy, sz), order= sOrder)
    save_high_res(highResB0Path, sp_high, lowResMaskHdr, b0HighRes)
    # save_high_res(highResB0Path, sp_high, lowResMaskHdr,
    #               applymask(b0HighRes, resize(np.double(b0_orig), (sx, sy, sz), order= sOrder)))

    # unring the b0
    check_call(['unring.a64', highResB0Path, highResB0Path])
    b0_gibs = nib.load(highResB0Path).get_data()

    ind_x, ind_y, ind_z, _= np.where(lowResImg>0)
    b0_gibs[b0_gibs > b0.max()] = b0[ind_x, ind_y, ind_z].max()
    b0_gibs[b0_gibs < b0.min()] = b0[ind_x, ind_y, ind_z].min()
    # save_high_res(highResB0Path, sp_high, lowResImgHdr, b0_gibs)


    # insert b0 back ------------------------------------------------------------------
    for i in where_b0:
        highResImg[:,:,:,i]= b0_gibs

    highResImg[highResImg > (lowResImg+eps).max()] = lowResImg[ind_x, ind_y, ind_z,: ].max()
    highResImg[highResImg < (lowResImg-eps).min()] = lowResImg[ind_x, ind_y, ind_z,: ].min()

    highResImg= applymask(highResImg, highResMask)
    highResImgPath= lowResImgPath.split('.')[0]+'_resampled'+'.nii.gz'
    save_high_res(highResImgPath, sp_high, lowResImgHdr, highResImg)

    return (highResImgPath, highResMaskPath)


if __name__=='__main__':
    from dipy.io import read_bvals_bvecs

    # lowResImgPath= '/home/tb571/Downloads/Harmonization-Python/connectom_prisma_demoData/A/connectom/dwi_A_connectom_st_b1200.nii.gz'
    # lowResMaskPath= '/home/tb571/Downloads/Harmonization-Python/connectom_prisma_demoData/A/connectom/mask.nii.gz'
    # bvals, _= read_bvals_bvecs(
    #     '/home/tb571/Downloads/Harmonization-Python/connectom_prisma_demoData/A/connectom/dwi_A_connectom_st_b1200.bval',
    #     None)

    lowResImgPath= '/home/tb571/Downloads/Harmonization-Python/test_data/test_a/connectom/connectom_a_dwi.nii.gz'
    lowResMaskPath= '/home/tb571/Downloads/Harmonization-Python/test_data/test_a/connectom/connectom_a_mask.nii.gz'
    bvals, _= read_bvals_bvecs(
        '/home/tb571/Downloads/Harmonization-Python/test_data/test_a/connectom/connectom_a_dwi.bval', None)

    lowRes = nib.load(lowResImgPath)
    lowResImg= lowRes.get_data()
    lowResImgHdr= lowRes.header

    lowRes = nib.load(lowResMaskPath)
    lowResMask= lowRes.get_data()
    lowResMaskHdr= lowRes.header

    resampling(lowResImgPath, lowResMaskPath, lowResImg, lowResImgHdr, lowResMask, lowResMaskHdr,
               np.array([1.5,1.5,1.5]), bvals)


