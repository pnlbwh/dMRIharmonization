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

import numpy as np
from conversion import read_bvals, read_imgs_masks
from findBshells import findBShells
from util import dirname, basename, pjoin, B0_THRESH

def determineNshm(bvalFile):

    print(f'Determining maximum possible order of spherical harmonics for {bvalFile}\n')

    bvals= np.array(read_bvals(bvalFile))
    N_b= len(np.where(bvals>B0_THRESH)[0])

    if N_b<6:
        raise ValueError(f'At least 6 gradients are necessary for each b-shell, b-shell has only {N_b}')
    elif N_b>=6 and N_b<15:
        N_shm = 2
    elif N_b>=15 and N_b<28:
        N_shm = 4
    elif N_b>=28 and N_b<45:
        N_shm = 6
    else:
        N_shm = 8

    print(f'Maximum possible order is {N_shm} for {N_b} non-zero gradients\n')

    return (N_shm, N_b)

def verifySingleShellNess(bvalFile):

    print(f'Verifying {bvalFile} is single shell ...\n')

    quantized_bvals= findBShells(bvalFile)

    if len(quantized_bvals)>2:
        raise ValueError(f'{bvalFile} is not single shell. Use https://github.com/pnlbwh/multi-shell-dMRIharmonization')


def verifyNshmForAll(csvFile, N_shm):

    for imgPath in read_imgs_masks(csvFile)[0]:
        directory = dirname(imgPath)
        prefix = basename(imgPath).split('.nii')[0]
        bvalFile = pjoin(directory, prefix + '.bval')
        verifySingleShellNess(bvalFile)
        verifyNshm(N_shm, bvalFile)


def verifyNshm(nshm, bvalFile):

    if nshm>=2:
        print(f'Verifying suitability of spherical harmonics order {nshm} for {bvalFile}\n')

    N_shm, N_b= determineNshm(bvalFile)
    if nshm>=2 and nshm>N_shm:
        raise ValueError(f'Order of spherical harmonics {nshm} is higher than possible with {N_b} gradients for {bvalFile}. '
                         'See README.md and reduce --nshm')
    elif nshm>=2:
        print(f'Spherical harmonics order {nshm} is suitable, continuing\n')



if __name__=='__main__':
    pass


