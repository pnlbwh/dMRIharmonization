#!/usr/bin/env python

import argparse
from plumbum.cmd import fslroi
from plumbum import FG
def main():

    parser = argparse.ArgumentParser(description='Write pruned nifti data')
    parser.add_argument('--nifti', type=str, required=True,
                        help='nifti file')

    parser.add_argument('--bval', type=str, required=True,
                    help='bval file')

    parser.add_argument('--bvec', type=str, required=True,
                    help='bvec file')

    parser.add_argument('--prefix', type=str, required=True,
                    help='output prefix')


    ind = [1,3,7,11,15,17,19,21,23,27,33]

    args = parser.parse_args()
    nifitFile= args.nifti
    bvalFile= args.bval
    bvecFile= args.bvec
    prefix= args.prefix

    bvals= read_bvals(bvalFile)
    bvecs= read_bvecs(bvecFile)

    bvals= [bvals[i] for i in ind]
    bvecs= [bvecs[i] for i in ind]

    write_bvals(prefix+'.bval', bvals)
    write_bvecs(prefix+'.bvec', bvecs)

    xmin = str(10)
    xsize = str(70)
    ymin = str(10)
    ysize = str(70)
    zmin = str(10)
    zsize = str(40)
    tmin = (',').join(str(i) for i in ind)
    tsize = str(len(ind))
    fslroi[nifitFile, prefix+'.nii.gz', xmin, xsize, ymin, ysize, zmin, zsize, tmin, tsize] & FG

def read_bvecs(bvec_file):

    with open(bvec_file, 'r') as f:
        bvecs = [[float(num) for num in line.split()] for line in f.read().split('\n') if line]

    # bvec_file can be 3xN or Nx3
    # we want to return as Nx3
    if len(bvecs) == 3:
        bvecs = tranpose(bvecs)

    return bvecs


def read_bvals(bval_file):

    with open(bval_file, 'r') as f:
        bvals = [float(num) for num in f.read().split()]

    # bval_file can be 1 line or N lines
    return bvals


def write_bvals(bval_file, bvals):
    with open(bval_file, 'w') as f:
        f.write(('\n').join(str(b) for b in bvals))


def write_bvecs(bvec_file, bvecs):
    with open(bvec_file, 'w') as f:
        # when bvecs is a list
        f.write(('\n').join((' ').join(str(i) for i in row) for row in bvecs))


def tranpose(bvecs):
    bvecs_T = list(map(list, zip(*bvecs)))

    return bvecs_T

if __name__=='__main__':
    main()

'''
A

/home/tb571/Downloads/Harmonization-Python/lib/test_data_creation.py \
--prefix /home/tb571/Downloads/Harmonization-Python/test_data/test_a/connectom/connectom_a_dwi \
--nifti /home/tb571/Downloads/Harmonization-Python/connectom_prisma_demoData/A/connectom/dwi_A_connectom_st_b1200.nii.gz \
--bval /home/tb571/Downloads/Harmonization-Python/connectom_prisma_demoData/A/connectom/dwi_A_connectom_st_b1200.bval \
--bvec /home/tb571/Downloads/Harmonization-Python/connectom_prisma_demoData/A/connectom/dwi_A_connectom_st_b1200.bvec

bet /home/tb571/Downloads/Harmonization-Python/test_data/test_a/connectom/connectom_a_dwi.nii.gz \
/home/tb571/Downloads/Harmonization-Python/test_data/test_a/connectom/connectom_a_mask -m

'''