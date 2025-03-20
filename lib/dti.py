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

def dti(imgPath, maskPath, inPrefix, outPrefix, tool='FSL'):

    vol = load(imgPath)
    mask = load(maskPath)
    masked_vol = applymask(vol.get_fdata(), mask.get_fdata())

    if tool=='DIPY':
        print('dipy dtifit ', imgPath)

        bvals, bvecs = read_bvals_bvecs(inPrefix + '.bval', inPrefix + '.bvec')

        gtab= gradient_table(bvals, bvecs, b0_threshold= B0_THRESH)
        dtimodel= dipyDti.TensorModel(gtab, fit_method="LS")
        dtifit= dtimodel.fit(masked_vol)
        fa= dtifit.fa
        md= dtifit.md

        save_nifti(outPrefix + '_FA.nii.gz', fa, vol.affine, vol.header)
        save_nifti(outPrefix + '_MD.nii.gz', md, vol.affine, vol.header)


    elif tool=='FSL':
        from plumbum.cmd import dtifit
        from plumbum import FG

        print('fsl dtifit ', imgPath)
        dtifit['-k', imgPath,
               '-m', maskPath,
               '-r', inPrefix + '.bvec',
               '-b', inPrefix + '.bval',
               '-o', outPrefix
              ] & FG


    rms = np.sqrt(np.mean(np.power(masked_vol,2), axis=3))
    ga = np.std(masked_vol, axis=3)/(rms+10e-9)
    ga3D = ga.reshape(masked_vol.shape[0], masked_vol.shape[1], masked_vol.shape[2])
    gfa_vol = np.nan_to_num(ga3D)

    save_nifti(outPrefix + '_GFA.nii.gz', gfa_vol, vol.affine, vol.header)


if __name__ == '__main__':
    pass

