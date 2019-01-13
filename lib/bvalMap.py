#!/usr/bin/env python
import warnings
with warnings.catch_warnings():
    warnings.filterwarnings("ignore", category=FutureWarning)

    from dipy.io.image import load_nifti, save_nifti
    from dipy.io import read_bvals_bvecs
    from dipy.reconst.shm import normalize_data
    from dipy.segment.mask import applymask
    from scipy.io import savemat

import numpy as np

def bvalMap(dwi, bvals, bmax):

    # find b0
    where_b0= np.where(bvals == 0)[0]
    b0= dwi[...,where_b0].mean(-1)

    # normalize dwi by b0
    dwiPrime= normalize_data(dwi, where_b0)

    # scale signal to the power of bmax/b
    ratio= []
    for b in bvals:
        if b:
            ratio.append(bmax/b)
        else:
            ratio.append(1.)

    ratio= np.reshape(ratio, (1, len(bvals)))
    dwiHat= dwiPrime**ratio

    # un-normalize dwi by b0
    dwiNew= applymask(dwiHat, b0)
    bvalsNew= bvals*bmax/max(bvals)

    return (dwiNew, bvalsNew)


if __name__=='__main__':
    dwi, affine= load_nifti('/home/tb571/Downloads/Harmonization-Python/connectom_prisma_demoData/A/connectom/dwi_A_connectom_st_b1200.nii.gz')
    mask= load_nifti('/home/tb571/Downloads/Harmonization-Python/connectom_prisma_demoData/A/connectom/mask.nii.gz')[0]
    bvals, _= read_bvals_bvecs(
        '/home/tb571/Downloads/Harmonization-Python/connectom_prisma_demoData/A/connectom/dwi_A_connectom_st_b1200.bval',
        '/home/tb571/Downloads/Harmonization-Python/connectom_prisma_demoData/A/connectom/dwi_A_connectom_st_b1200.bvec')
    dwi= bvalMap(dwi, bvals, 300)
    savemat('/home/tb571/Downloads/Harmonization-Python/connectom_prisma_demoData/bmap/b_data.mat',{'dwi_mat':dwi})
    save_nifti('/home/tb571/Downloads/Harmonization-Python/connectom_prisma_demoData/bmap/dwi.nii.gz', dwi, affine)

    '''
    # for testing in MATLAB
        
    dwPrime= np.random.rand(3*3*5*5).reshape((3,3,5,5))
    bvals= np.reshape([50, 100, 200, 50, 200], (5,1))
    bmax= 300
    ratio= []
    for b in bvals:
        ratio.append(bmax/b)

    dwiHat = dwiPrime ** ratio
    from scipy.io import savemat
    savemat('test_dwi.mat',{'S':dwiPrime, 'dwiHat':dwiHat, 'b_o':bvals, 'b_n':bmax})
    
    '''