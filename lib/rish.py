#!/usr/bin/env python

import warnings
with warnings.catch_warnings():
    warnings.filterwarnings("ignore", category=FutureWarning)

    from dipy.io.image import load_nifti, save_nifti
    from dipy.io import read_bvals_bvecs
    from dipy.core.gradients import gradient_table
    from dipy.reconst.shm import QballModel
    from dipy.segment.mask import applymask

from normalize import normalize_data, find_b0

import numpy as np
import os

def rish(imgPath, maskPath, inPrefix, outPrefix, N_shm, qb_model= None):

    data, affine= load_nifti(imgPath)
    mask_data, _= load_nifti(maskPath)

    if not qb_model:
        print('Computing shm_coeff of ', imgPath)
        bvals, bvecs = read_bvals_bvecs(inPrefix+'.bval', inPrefix+'.bvec')

        gtab = gradient_table(bvals,  bvecs)
        qb_model = QballModel(gtab, sh_order=N_shm)

        # save baseline image
        b0 = find_b0(data, where_b0=np.where(qb_model.gtab.b0s_mask)[0])
        if not os.path.exists(inPrefix+'_bse.nii.gz'):
            save_nifti(inPrefix+'_bse.nii.gz', applymask(b0, mask_data), affine= affine)
    else:
        b0= None


    # inserting correct shm_coeff computation block ---------------------------------
    smooth= 0.006
    data = applymask(data, mask_data)
    data_norm, _ = normalize_data(data, where_b0=np.where(qb_model.gtab.b0s_mask)[0])


    L= qb_model.n*(qb_model.n+1)
    L**=2
    _fit_matrix= np.linalg.pinv(qb_model.B.T @ qb_model.B+ np.diag(smooth*L)) @ qb_model.B.T
    shm_coeff= np.dot(data_norm[..., qb_model._where_dwi], _fit_matrix.T)
    shm_coeff= applymask(shm_coeff, mask_data)
    # -------------------------------------------------------------------------------

    shm_coeff_squared= shm_coeff**2
    shs_same_level= [[0, 1], [1, 6], [6, 15], [15, 28], [28, 45]]
    # rishImgs= np.zeros((data.shape[0], data.shape[1], data.shape[2], N_shm), dtype= float)
    for i in range(0, N_shm+1, 2):
        ind= int(i/2)
        temp= np.sum(shm_coeff_squared[:,:,:,shs_same_level[ind][0]:shs_same_level[ind][1]], axis= 3)
        save_nifti(f'{outPrefix}_L{ind*2}.nii.gz', temp, affine)

        # rishImgs[:, :, :, i]= temp


    return (b0, shm_coeff, qb_model)

if __name__=='__main__':
    pass