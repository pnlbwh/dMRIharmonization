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

def find_b0(dwi, where_b0, mask= None):
    b0= dwi[...,where_b0].mean(-1)
    np.nan_to_num(b0).clip(min= 0., out= b0)

    if mask is not None:
        return b0*mask
    else:
        return b0


def normalize_data(dwi, where_b0= None, mask= None, b0= None):

    dwi= dwi.astype('float32')

    if where_b0 is not None and b0 is None:
        b0= find_b0(dwi, where_b0, mask)
        np.nan_to_num(b0).clip(min=1., out=b0)
        for i in where_b0:
            dwi[...,i] = b0
    else:
        np.nan_to_num(b0).clip(min=1., out=b0)


    dwiPrime= dwi/b0[...,None]
    np.nan_to_num(dwiPrime).clip(min=0., max=1., out= dwiPrime)

    if mask is not None:
        dwiPrime= applymask(dwiPrime, mask)

    return (dwiPrime, b0)