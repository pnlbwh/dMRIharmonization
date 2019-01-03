import warnings

with warnings.catch_warnings():
    warnings.filterwarnings("ignore", category=FutureWarning)

    from dipy.io.image import load_nifti, save_nifti
    from dipy.io import read_bvals_bvecs
    from dipy.core.gradients import gradient_table
    from dipy.reconst.shm import QballModel
    from dipy.segment.mask import applymask

import numpy as np
eps= 2.2204e-16

def rish(imgPath, maskPath, inPrefix, outPrefix, N_shm, qb_model= None):

    data, affine= load_nifti(imgPath)
    mask_data, _= load_nifti(maskPath)

    if not qb_model:
        bvals, bvecs = read_bvals_bvecs(inPrefix+'.bval', inPrefix+'.bvec')

        gtab = gradient_table(bvals,  bvecs)
        qb_model = QballModel(gtab, sh_order=N_shm)

        # save baseline image
        b0 = data[..., np.where(bvals==0)[0]].mean(-1)
        save_nifti(inPrefix+'_bse.nii.gz', applymask(b0, mask_data), affine= affine)
    else:
        b0= None

    qb_fit = qb_model.fit(data, mask_data)
    shm_coeff= qb_fit.shm_coeff

    shm_coeff_squared= shm_coeff**2
    shs_same_level= [[0, 1], [1, 6], [6, 15], [15, 28], [28, 45]]
    # rishImgs= np.zeros((data.shape[0], data.shape[1], data.shape[2], N_shm), dtype= float)
    for i in range(0, N_shm+1, 2):
        ind= int(i/2)
        temp= np.sum(shm_coeff_squared[:,:,:,shs_same_level[ind][0]:shs_same_level[ind][1]], axis= 3)
        save_nifti(f'{outPrefix}_L{ind*2}.nii.gz', temp, affine)

        # rishImgs[:, :, :, i]= temp

    return (b0, shm_coeff, qb_model)

if __name__ == '__main__':
    rish('/home/tb571/Downloads/Harmonization-Python/connectom_prisma_demoData/A/connectom/dwi_A_connectom_st_b1200.nii.gz',
        '/home/tb571/Downloads/Harmonization-Python/connectom_prisma_demoData/A/connectom/mask.nii.gz',
        '/home/tb571/Downloads/Harmonization-Python/connectom_prisma_demoData/A/connectom/dwi_A_connectom_st_b1200',
        '/home/tb571/Downloads/Harmonization-Python/connectom_prisma_demoData/py_thon/dwi',
         8)