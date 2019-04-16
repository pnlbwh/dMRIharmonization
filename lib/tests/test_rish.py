#!/usr/bin/env python

from test_util import *
from rish import rish

class TestRish(unittest.TestCase):

    def test_rish_feature(self):
        inPath= pjoin(FILEDIR, 'connectom_prisma/connectom/A/')
        outPath= pjoin(inPath, 'harm')

        if exists(outPath):
            rmtree(outPath)
        mkdir(outPath)
        inPrefix= pjoin(inPath, 'dwi_A_connectom_st_b1200')
        outPrefix= pjoin(outPath, 'dwi_A_connectom_st_b1200')

        dwiPath= inPrefix+'.nii.gz'
        maskPath= inPrefix+'_mask.nii.gz'
        rish(dwiPath, maskPath, inPrefix, outPrefix, 6)


if __name__ == '__main__':
    unittest.main()


