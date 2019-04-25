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

import os, sys
import warnings
with warnings.catch_warnings():
    warnings.filterwarnings("ignore", category=FutureWarning)
    import nibabel as nib
    from dipy.io.image import save_nifti
import numpy as np

np.set_printoptions(precision=5)

def main():

    if sys.argv[1]=='-h' or sys.argv[1]=='--help':
        print(
'''
This module is the gateway for testing equivalence between MATALAB and PYTHON rish features.
Example usage:
./rish_diff_map.py abs/path/to/py_rish/prefix abs/path/to/mat_rish/prefix 6

rish features are loaded from path/to/py_rish_L{i} for i=0:2:N_shm, likewise for mat_rish
Specify a sample path in 'imgPath_py' and 'imgPath_mat' using *XYZ* for case, *ORDER* for shm order in sample path.
The program will replace *XYZ* with caseid obtained from provided caselist. It will also substitute *ORDER* with
proper spherical harmonic order upto N_shm given to
rish_diff(imgPath_given_mat, imgPath_given_py, caselist, N_shm= 6) and
scale_diff(imgPath_given_mat, imgPath_given_py, N_shm= 6)
'''
        )
        exit()


    N_shm=int(sys.argv[3])

    dir_bak= os.getcwd()

    directory= os.path.dirname(sys.argv[1])
    prefix= os.path.basename(sys.argv[1])
    os.chdir(directory)
    py=[]
    for i in range(0,N_shm+1,2):
        img= nib.load(f'{prefix}_L{i}.nii.gz')
        py.append(img.get_data())

    pyX, pyY, pyZ = np.shape(py)[1: ]


    directory= os.path.dirname(sys.argv[2])
    prefix= os.path.basename(sys.argv[2])
    os.chdir(directory)
    mat=[]
    for i in range(0,N_shm+1,2):
        img= nib.load(f'{prefix}_L{i}.nii.gz')
        mat.append(img.get_data()[:pyX, :pyY, :pyZ])

        ind= int(i/2)
        save_nifti(f'diff_L{ind}.nii.gz', py[ind]-mat[ind], affine= img.affine, header= img.header)


    mat= np.array(mat)
    py= np.array(py)

    for i in range(0,int(N_shm/2)+1):

        temp_py= py[i].flatten()
        temp_mat= mat[i].flatten()
        print(f'Rish L{i}\n')
        print(f'# of non zero voxels: '
              f'python {len(np.where(temp_py!=0)[0])}\t\t'
              f'matlab {len(np.where(temp_mat!=0)[0])}')

        print(f'Maximum values: python {temp_py.max()}\t matlab {temp_mat.max()}')
        print(f'Minimum values: python {temp_py.min()}\t\t\t matlab {temp_mat.min()}')

        diff= abs(temp_py- temp_mat)
        print('\n')
        print(f'Mean difference over all the voxels: {diff.mean()}')
        print(f'Standard deviation over all the voxels: {diff.std()}')

        print('\n\n')

    os.chdir(dir_bak)

if __name__=='__main__':
    main()
