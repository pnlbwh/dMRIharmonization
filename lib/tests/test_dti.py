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
from dti import dti


class TestDti(unittest.TestCase):

    def test_dti(self):
        inPath = pjoin(FILEDIR, 'connectom_prisma/connectom/A/')
        outPath = pjoin(inPath, 'dti')

        if exists(outPath):
            rmtree(outPath)
        mkdir(outPath)
        inPrefix = pjoin(inPath, 'dwi_A_connectom_st_b1200')

        lowResImgPath = inPrefix + '.nii.gz'
        lowResMaskPath = inPrefix + '_mask.nii.gz'

        outPrefix = pjoin(outPath, 'dwi_A_connectom_st_b1200')

        dti(lowResImgPath, lowResMaskPath, inPrefix, outPrefix)


if __name__ == '__main__':
    unittest.main()
