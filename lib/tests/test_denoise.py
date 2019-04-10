#!/usr/bin/env python

from util import *
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
        save_nifti(outPrefix + '.nii.gz', dwiNew, dwi.affine)
        copyfile(inPrefix + '.bvec', outPrefix + '.bvec')
        copyfile(inPrefix + '.bval', outPrefix + '.bval')


if __name__ == '__main__':
    unittest.main()