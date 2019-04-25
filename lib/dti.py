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

from util import *

def dti(imgPath, maskPath, inPrefix, outPrefix):

    print('dipy dtifit ', imgPath)

    vol= load(imgPath)
    mask= load(maskPath)
    bvals, bvecs= read_bvals_bvecs(inPrefix + '.bval', inPrefix + '.bvec')
    masked_vol = applymask(vol.get_data(), mask.get_data())

    gtab= gradient_table(bvals, bvecs)
    dtimodel= dipyDti.TensorModel(gtab, fit_method="LS")
    dtifit= dtimodel.fit(masked_vol)
    fa= dtifit.fa
    md= dtifit.md

    save_nifti(outPrefix + '_FA.nii.gz', fa, vol.affine, vol.header)
    save_nifti(outPrefix + '_MD.nii.gz', md, vol.affine, vol.header)

    gfa_vol= gfa(masked_vol)
    save_nifti(outPrefix+'_GFA.nii.gz', gfa_vol, vol.affine, vol.header)


if __name__ == '__main__':
    pass