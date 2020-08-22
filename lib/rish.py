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

from normalize import normalize_data, find_b0
from util import *


def rish(imgPath, maskPath, inPrefix, outPrefix, N_shm, qb_model=None):

    img = load(imgPath)
    data = img.get_data()
    affine = img.affine
    hdr = img.header
    mask_data = load(maskPath).get_data()

    if not qb_model:
        print('Computing shm_coeff of ', imgPath)
        bvals, bvecs = read_bvals_bvecs(inPrefix + '.bval', inPrefix + '.bvec')

        # make bvals and bvecs full sampled to encounter reconstruction error
        # in reconstSignal.py
        bvals = np.append(bvals, bvals)
        bvecs = np.append(bvecs, -bvecs, axis=0)
        data = np.append(data, data, axis=3)

        gtab = gradient_table(bvals, bvecs, b0_threshold=B0_THRESH)
        qb_model = QballModel(gtab, sh_order=N_shm)

        # save baseline image
        b0 = find_b0(data, where_b0=np.where(qb_model.gtab.b0s_mask)[0])
        if not isfile(inPrefix + '_bse.nii.gz'):
            save_nifti(
                inPrefix +
                '_bse.nii.gz',
                applymask(
                    b0,
                    mask_data),
                affine,
                hdr)
    else:
        b0 = None

    # inserting correct shm_coeff computation block --------------------------
    smooth = 0.00001
    data = applymask(data, mask_data)
    data_norm, _ = normalize_data(
        data, where_b0=np.where(
            qb_model.gtab.b0s_mask)[0])

    L = qb_model.n * (qb_model.n + 1)
    L **= 2
    _fit_matrix = np.linalg.pinv(
        qb_model.B.T @ qb_model.B +
        np.diag(
            smooth *
            L)) @ qb_model.B.T
    shm_coeff = np.dot(data_norm[..., qb_model._where_dwi], _fit_matrix.T)
    shm_coeff = applymask(shm_coeff, mask_data)
    # -------------------------------------------------------------------------------

    shm_coeff_squared = shm_coeff**2
    shs_same_level = [[0, 1], [1, 6], [6, 15], [15, 28], [28, 45]]
    for i in range(0, N_shm + 1, 2):
        ind = int(i / 2)
        temp = np.sum(
            shm_coeff_squared[:, :, :, shs_same_level[ind][0]:shs_same_level[ind][1]], axis=3)
        save_nifti(f'{outPrefix}_L{ind*2}.nii.gz', temp, affine, hdr)

    return (b0, shm_coeff, qb_model)


if __name__ == '__main__':
    pass
