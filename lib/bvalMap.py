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

from util import *
from normalize import normalize_data

def remapBval(dwi, mask, bvals, bNew):

    # find b0
    where_b0= np.where(bvals <= B0_THRESH)[0]

    # normalize dwi by b0
    dwiPrime, b0 = normalize_data(dwi, where_b0, mask)
    # dwiPrime is already masked

    # scale signal to the power of bNew/b
    ratio= []
    bvalsNew= []
    for b in bvals:
        if b > B0_THRESH:
            bvalsNew.append(bNew)
            if (b>1.01*bNew or b<0.99*bNew):
                ratio.append(bNew / b)
            else:
                ratio.append(1.)
        else:
            ratio.append(1.)
            bvalsNew.append(0)

    ratio= np.reshape(ratio, (1, len(bvals)))
    dwiHat= dwiPrime**ratio

    # un-normalize dwi by b0
    dwiNew= applymask(dwiHat, b0)

    return (dwiNew, np.array(bvalsNew))


if __name__=='__main__':
    pass
