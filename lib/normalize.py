import numpy as np

def find_b0(dwi, where_b0):
    b0= dwi[...,where_b0].mean(-1)
    np.nan_to_num(b0).clip(min= 1., out= b0)

    return b0


def normalize_data(dwi, where_b0= None, b0= None):
    if b0 is None:
        b0= find_b0(dwi, where_b0)

    dwiPrime= dwi/b0[...,None]
    np.nan_to_num(dwiPrime).clip(min=0., max= 1., out= dwiPrime)
    # dwiPrime[np.isnan(dwiPrime)]= 0.
    # dwiPrime[dwiPrime<0]= 0.
    # dwiPrime[dwiPrime>1]= 1.

    return (dwiPrime, b0)