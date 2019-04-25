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

from test_util import *
from denoising import denoising
from preprocess import nrrd2nifti


class TestDenoise(unittest.TestCase):

    def test_denoise(self):
        inPath = pjoin(FILEDIR, 'connectom_prisma/connectom/A/')
        inPrefix = pjoin(inPath, 'dwi_A_connectom_st_b1200')

        lowResImgPath = inPrefix + '.nii.gz'
        lowResMaskPath = inPrefix + '_mask.nii.gz'

        # load signal attributes for pre-processing ----------------------------------------------------------------
        imgPath = nrrd2nifti(lowResImgPath)
        dwi = load(imgPath)

        maskPath = nrrd2nifti(lowResMaskPath)
        mask = load(maskPath)

        print('Denoising ', imgPath)
        dwiNew, _= denoising(dwi.get_data(), mask.get_data())
        outPrefix = imgPath.split('.')[0] + '_denoised'
        save_nifti(outPrefix + '.nii.gz', dwiNew, dwi.affine, dwi.header)
        copyfile(inPrefix + '.bvec', outPrefix + '.bvec')
        copyfile(inPrefix + '.bval', outPrefix + '.bval')


if __name__ == '__main__':
    unittest.main()