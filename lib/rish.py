import warnings

with warnings.catch_warnings():
    warnings.filterwarnings("ignore", category=FutureWarning)

    from dipy.io.image import load_nifti, save_nifti
    from dipy.io import read_bvals_bvecs
    from dipy.core.gradients import gradient_table
    from dipy.reconst.shm import QballModel

import numpy as np

def rish(imgPath, maskPath, prefix, outPrefix):

    data, affine= load_nifti(imgPath)
    mask_data, _= load_nifti(maskPath)
    bvals, bvecs = read_bvals_bvecs(prefix+'.bval', prefix+'.bvec')

    gtab = gradient_table(bvals,  bvecs)
    shm_order= 8
    qb_model = QballModel(gtab, sh_order=shm_order)
    qb_fit = qb_model.fit(data, mask_data)
    shm_coeff= qb_fit.shm_coeff
    shm_coeff= shm_coeff**2

    shs_same_level= [[0, 1], [1, 6], [6, 15], [15, 28], [28, 45]]
    N_shm= int(shm_order/2)+1

    # rishImgs= np.zeros((data.shape[0], data.shape[1], data.shape[2], N_shm), dtype= float)
    for i in range(N_shm):

        temp= np.sum(shm_coeff[:,:,:,shs_same_level[i][0]:shs_same_level[i][1]], axis= 3)
        save_nifti(f'{outPrefix}_{i*2}.nii.gz', temp, affine, None)

        # rishImgs[:, :, :, i]= temp