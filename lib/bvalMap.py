#!/usr/bin/env python
import warnings
with warnings.catch_warnings():
    warnings.filterwarnings("ignore", category=FutureWarning)

    from dipy.io.image import load_nifti, save_nifti
    from dipy.io import read_bvals_bvecs
    from normalize import normalize_data
    from dipy.segment.mask import applymask
    from scipy.io import savemat

import numpy as np

def remapBval(dwi, bvals, bNew):

    # find b0
    where_b0= np.where(bvals == 0)[0]

    # normalize dwi by b0
    dwiPrime, b0= normalize_data(dwi, where_b0)

    # scale signal to the power of bNew/b
    ratio= []
    bvalsNew= []
    for b in bvals:
        if b:
            bvalsNew.append(bNew)
            if (b>1.01*bNew or b<0.99*bNew):
                ratio.append(bNew / b)
            else:
                ratio.append(1.)
        else:
            ratio.append(1.)
            bvalsNew.append(0)

    ratio= np.reshape(ratio, (1, len(bvals)))
    dwiHat= dwiPrime**ratio

    # un-normalize dwi by b0
    dwiNew= applymask(dwiHat, b0)

    return (dwiNew, np.array(bvalsNew))


if __name__=='__main__':
    dwi, affine= load_nifti('/home/tb571/Downloads/Harmonization-Python/connectom_prisma_demoData/A/connectom/dwi_A_connectom_st_b1200.nii.gz')
    mask= load_nifti('/home/tb571/Downloads/Harmonization-Python/connectom_prisma_demoData/A/connectom/mask.nii.gz')[0]
    bvals, _= read_bvals_bvecs(
        '/home/tb571/Downloads/Harmonization-Python/connectom_prisma_demoData/A/connectom/dwi_A_connectom_st_b1200.bval',
        None)
    dwi, bvalsNew= remapBval(dwi, bvals, 300)
    # save_nifti('/home/tb571/Downloads/Harmonization-Python/connectom_prisma_demoData/bmap/dwi.nii.gz', dwi, affine)

    '''
    # for testing in MATLAB
    # savemat('/home/tb571/Downloads/Harmonization-Python/connectom_prisma_demoData/bmap/b_data.mat',{'dwi_mat':dwi})
    dwPrime= np.random.rand(3*3*5*5).reshape((3,3,5,5))
    bvals= np.reshape([50, 100, 200, 50, 200], (5,1))
    bNew= 300
    ratio= []
    for b in bvals:
        ratio.append(bNew/b)

    dwiHat = dwiPrime ** ratio
    from scipy.io import savemat
    savemat('test_dwi.mat',{'S':dwiPrime, 'dwiHat':dwiHat, 'b_o':bvals, 'b_n':bNew})
    
    '''