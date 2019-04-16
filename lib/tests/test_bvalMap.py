#!/usr/bin/env python

from test_util import *
from bvalMap import remapBval
from preprocess import nrrd2nifti
from dipy.io import read_bvals_bvecs

class TestBmap(unittest.TestCase):

    def test_bvalMap(self):
        inPath= pjoin(FILEDIR, 'connectom_prisma/connectom/A/')
        inPrefix= pjoin(inPath, 'dwi_A_connectom_st_b1200')

        lowResImgPath= inPrefix+'.nii.gz'
        bvalPath= inPrefix+'.bval'
        lowResMaskPath= inPrefix+'_mask.nii.gz'

        # load signal attributes for pre-processing ----------------------------------------------------------------
        imgPath = nrrd2nifti(lowResImgPath)
        dwi = load(imgPath)

        maskPath = nrrd2nifti(lowResMaskPath)
        mask = load(maskPath)

        bvals, _ = read_bvals_bvecs(bvalPath, None)

        bNew= 1000.

        print('B value mapping ', imgPath)
        dwiNew, bvalsNew= remapBval(dwi.get_data(), mask.get_data(), bvals, bNew)

        outPrefix = imgPath.split('.')[0] + '_bmapped'
        save_nifti(outPrefix + '.nii.gz', dwiNew, dwi.affine, dwi.header)
        copyfile(inPrefix + '.bvec', outPrefix + '.bvec')
        write_bvals(outPrefix + '.bval', bvals)

if __name__ == '__main__':
    unittest.main()