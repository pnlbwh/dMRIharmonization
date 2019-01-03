#!/usr/bin/env python

import os, sys
import warnings
with warnings.catch_warnings():
    warnings.filterwarnings("ignore", category=FutureWarning)
    import nibabel as nib
    from dipy.io.image import save_nifti
import numpy as np

np.set_printoptions(precision=5)

def main():

    py=sys.argv[1]
    mat=sys.argv[2]
    N_shm=int(sys.argv[3])

    dir_bak= os.getcwd()

    os.chdir(sys.argv[1])
    py=[]
    for i in range(0,N_shm+1,2):
        py.append(nib.load(f'dwi_L{i}.nii.gz').get_data())

    os.chdir(sys.argv[2])
    mat=[]
    for i in range(0,N_shm+1,2):
        img= nib.load(f'dwi_L{i}.nii.gz')
        mat.append(img.get_data())

        ind= int(i/2)
        save_nifti(f'diff_L{ind}.nii.gz', py[ind]-mat[ind], affine= img.affine)


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
