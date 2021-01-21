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

from conversion import read_bvals, write_bvals
import numpy as np
from os.path import abspath
from util import B0_THRESH, B_QUANT, BSHELL_MIN_DIST


def usage():
    print('''Find b-shells in a DWI 
Usage:
findBShells /path/to/input/bval/file /output/file/to/write/bshells
''')


def findBShells(bvalFile, outputBshellFile= None):

    given_bvals= read_bvals(abspath(bvalFile))

    # get unique bvalues in ascending order
    unique_bvals= np.unique(given_bvals)

    # identify b0s
    quantized_bvals= unique_bvals.copy()
    quantized_bvals[unique_bvals<=B0_THRESH]= 0.

    # round to multiple of B_QUANT (50 or 100)
    quantized_bvals= np.unique(np.round(quantized_bvals/B_QUANT)*B_QUANT)

    print('b-shell bvalues', quantized_bvals)

    for bval in quantized_bvals:
        print('Indices corresponding to b-shell', bval)
        print(np.where(abs(bval-given_bvals)<=BSHELL_MIN_DIST)[0],'\n')


    if outputBshellFile:
        print('Saving the b-shell bvalues in', outputBshellFile)
        write_bvals(outputBshellFile, quantized_bvals)

    return quantized_bvals

if __name__== '__main__':
    import sys
    if len(sys.argv)==1 or sys.argv[1]=='-h' or sys.argv[1]=='--help':
        usage()
        exit()

    findBShells(sys.argv[1],sys.argv[2])

