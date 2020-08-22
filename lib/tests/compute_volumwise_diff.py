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

import warnings
with warnings.catch_warnings():
    warnings.filterwarnings("ignore", category=FutureWarning)
    import nibabel as nib

from conversion import read_bvals, read_bvecs
import numpy as np
import sys
from test_util import B0_THRESH


np.set_printoptions(precision=5)
# BINS= [0, 1, 5, 10, np.inf]
# BINS= [0, 1, 10, 50, 100]
# BINS= [0, 5, 10, 50, np.inf]
BINS = [0, 5, 10, 50, 100]


def hist_calc(a, bins):

    a = np.array([item for sublist in a for item in sublist])

    hist_string = []
    N = len(bins)
    for i in range(N - 1):
        hist_string.append(f'{bins[i]} <--> {bins[i+1]}')

    temp = np.histogram(a, bins=bins)[0]
    hist = temp.astype('float') / a.size

    try:
        hist[np.isnan] = 0
    except BaseException:
        pass

    print('%20s : %s' % ('Bins', 'Density'))
    for i in range(N - 1):
        print('%20s : %.6f' % (hist_string[i], hist[i]))

    return hist

# bval/bvec difference


def bval_bvec_difference(imgPath_given_mat, imgPath_given_py, caselist):
    bval_diff_list = []
    bval_total_diff = []
    bval_mean_diff = []

    bvec_diff_list = []
    bvec_total_diff = []
    bvec_mean_diff = []

    for case in caselist:
        imgPath_mat = imgPath_given_mat.replace('XYZ', case)
        imgPath_py = imgPath_given_py.replace('XYZ', case)

        # load python
        bvals_py = np.array(read_bvals(imgPath_py.split('.nii')[0] + '.bval'))
        bvecs_py = np.array(read_bvecs(imgPath_py.split('.nii')[0] + '.bvec'))

        bvals_py, bvecs_py = stack_one_b0(bvals_py, bvecs=bvecs_py)
        bvals_py = bvals_py.astype('float64')
        bvecs_py = bvecs_py.astype('float64')

        # load matlab
        bvals_mat = np.array(
            read_bvals(
                imgPath_mat.split('.nii')[0] +
                '.bval'))
        bvecs_mat = np.array(
            read_bvecs(
                imgPath_mat.split('.nii')[0] +
                '.bvec'))
        bvals_mat = bvals_mat.astype('float64')
        bvecs_mat = bvecs_mat.astype('float64')

        # obtain case wise difference for bvals
        diff = 2 * abs(bvals_py - bvals_mat) / (bvals_py + bvals_mat) * 100
        bval_diff_list.append(diff.flatten())
        bval_mean_diff.append(diff.mean())

        # obtain case wise difference
        diff = 2 * abs(bvecs_py - bvecs_mat) / (bvals_py + bvals_mat) * 100
        bvec_diff_list.append(diff.flatten())
        bvec_mean_diff.append(diff.mean())

    print(
        f'Mean of mean relative percentage difference over all the bvals: {np.mean(bval_mean_diff)}')
    hist_calc(bval_diff_list, bins=BINS)
    print('\n\n')

    print(
        f'Mean of mean relative percentage difference over all the bvecs: {np.mean(bvec_mean_diff)}')
    hist_calc(bvec_diff_list, bins=BINS)
    print('\n\n')


# resampled difference
# harmonized difference
# load python image, keep only the first b0
def dwi_difference(imgPath_given_mat, imgPath_given_py, caselist):
    diff_list = []
    total_diff = []
    mean_diff = []

    for case in caselist:
        imgPath_mat = imgPath_given_mat.replace('XYZ', case)
        imgPath_py = imgPath_given_py.replace('XYZ', case)

        # load python
        img_py_obj = nib.load(imgPath_py)
        img_py = img_py_obj.get_data()
        # keep only the first b0
        bvals = np.array(read_bvals(imgPath_py.split('.nii')[0] + '.bval'))
        img_py = stack_one_b0(bvals, dwi=img_py)
        img_py = img_py.astype('float64')

        pyX, pyY, pyZ, _ = img_py.shape

        # load matlab
        img_mat = nib.load(imgPath_mat).get_data()[:pyX, :pyY, :pyZ, :]
        img_mat = img_mat.astype('float64')

        # obtain case wise difference
        diff = 2 * abs(img_py - img_mat) / \
            (img_py + img_mat).clip(min=1.) * 100
        diff_list.append(diff.flatten())
        mean_diff.append(diff.mean())

        # save the difference mask
        nib.Nifti1Image(
            diff,
            img_py_obj.affine,
            img_py_obj.header).to_filename(
            imgPath_mat.split('.nii')[0] +
            '_diff_mask.nii.gz')

    print(
        f'Mean of mean relative percentage difference over all the voxels: {np.mean(mean_diff)}')
    hist_calc(diff_list, bins=BINS)


def stack_one_b0(bvals, bvecs=None, dwi=None):

    bvals = np.array(bvals)
    where_b0 = np.where(bvals <= B0_THRESH)[0]
    where_dwi = np.where(bvals > B0_THRESH)[0]

    reduced_ind = [where_b0[0]]
    for i in where_dwi:
        reduced_ind.append(i)

    if dwi is not None:
        reduced_volume = np.take(dwi, indices=reduced_ind, axis=3)
        return reduced_volume
    else:
        reduced_bvals = bvals[reduced_ind]
        reduced_bvecs = bvecs[reduced_ind, :]
        return (reduced_bvals, reduced_bvecs)


# dti difference
def dti_diff(imgPath_given_mat, imgPath_given_py, caselist):

    diff_list = []
    total_diff = []
    mean_diff = []

    for case in caselist:
        imgPath_mat = imgPath_given_mat.replace('XYZ', case)
        imgPath_py = imgPath_given_py.replace('XYZ', case)

        # load python
        img_py_obj = nib.load(imgPath_py)
        img_py = img_py_obj.get_data()
        img_py = img_py.astype('float64')
        pyX, pyY, pyZ = img_py.shape

        # load matlab
        img_mat = nib.load(imgPath_mat).get_data()[:pyX, :pyY, :pyZ]
        img_mat = img_mat.astype('float64')

        # obtain case wise difference
        diff = 2 * abs(img_py - img_mat) / \
            (img_py + img_mat).clip(min=1.) * 100
        diff_list.append(diff.flatten())
        mean_diff.append(diff.mean())

        # save the difference mask
        nib.Nifti1Image(
            diff,
            img_py_obj.affine,
            img_py_obj.header).to_filename(
            imgPath_mat.split('.nii')[0] +
            '_diff_mask.nii.gz')

    print(f'Mean of mean relative percentage difference: {np.mean(mean_diff)}')
    hist_calc(diff_list, bins=BINS)
    print('\n\n')

# rish difference


def rish_diff(imgPath_given_mat, imgPath_given_py, caselist, N_shm=6):

    for i in range(0, N_shm + 1, 2):

        imgPath_mat = imgPath_given_mat.replace('ORDER', str(i))
        imgPath_py = imgPath_given_py.replace('ORDER', str(i))

        diff_list = []
        total_diff = []
        mean_diff = []

        for case in caselist:
            imgPath_mat = imgPath_mat.replace('XYZ', case)
            imgPath_py = imgPath_py.replace('XYZ', case)

            # load python
            img_py_obj = nib.load(imgPath_py)
            img_py = img_py_obj.get_data()
            img_py = img_py.astype('float64')
            pyX, pyY, pyZ = img_py.shape

            # load matlab
            img_mat = nib.load(imgPath_mat).get_data()[:pyX, :pyY, :pyZ]
            img_mat = img_mat.astype('float64')

            # obtain case wise difference
            diff = 2 * abs(img_py - img_mat) / \
                (img_py + img_mat).clip(min=1.) * 100
            diff_list.append(diff.flatten())
            mean_diff.append(diff.mean())

            # save the difference mask
            affine = nib.load(imgPath_py).affine
            nib.Nifti1Image(
                diff,
                img_py_obj.affine,
                img_py_obj.header).to_filename(
                imgPath_mat.split('.nii')[0] +
                '_diff_mask.nii.gz')

        print(f'Feature: L{i}: ')
        print(
            f'Mean of mean relative percentage difference over all the voxels: {np.mean(mean_diff)}')
        hist_calc(diff_list, bins=BINS)
        print('\n\n')


# scale difference
def scale_diff(imgPath_given_mat, imgPath_given_py, N_shm=6):

    for i in range(0, N_shm + 1, 2):

        imgPath_mat = imgPath_given_mat.replace('ORDER', str(i))
        imgPath_py = imgPath_given_py.replace('ORDER', str(i))

        diff_list = []
        total_diff = []
        mean_diff = []

        # load python
        img_py_obj = nib.load(imgPath_py)
        img_py = img_py_obj.get_data()
        img_py = img_py.astype('float64')
        pyX, pyY, pyZ = img_py.shape

        # load matlab
        img_mat = nib.load(imgPath_mat).get_data()[:pyX, :pyY, :pyZ]
        img_mat = img_mat.astype('float64')

        # obtain case wise difference
        diff = 2 * abs(img_py - img_mat) / \
            (img_py + img_mat).clip(min=1.) * 100
        diff_list.append(diff.flatten())
        mean_diff.append(diff.mean())

        # save the difference mask
        affine = nib.load(imgPath_py).affine
        nib.Nifti1Image(
            diff,
            img_py_obj.affine,
            img_py_obj.header).to_filename(
            imgPath_mat.split('.nii')[0] +
            '_diff_mask.nii.gz')

        print(f'Feature: L{i}: ')
        print(
            f'Mean of mean relative percentage difference over all the voxels: {np.mean(mean_diff)}')
        hist_calc(diff_list, bins=BINS)
        print('\n\n')


def main():

    # default shm order
    N_shm = 4

    if len(sys.argv) == 2:
        if sys.argv[1] == '-h' or sys.argv[1] == '--help':
            print(
                '''
This module is the gateway for testing equivalence between MATALAB and PYTHON harmonization results.
Edit 'compute_volumwise_diff.py' for computing volumewise difference over all cases in a site.
Specify a sample path in 'imgPath_py' and 'imgPath_mat'. Use XYZ for caseid, ORDER for shm order in the sample paths.
The program will replace XYZ with caseids obtained from provided caselist. It will also substitute ORDER with
even shm order <= N_shm which is used to compute Rish and Scale differences for each shm.
Sample usage:
./compute_volumewise_diff.py -h   # print this usage
./compute_volumewise_diff.py 6    # N_shm=6 for Scale/Rish differences
./compute_volumeise_diff.py       # no need to specify shm order for DWI/Template/DTI differences
'''
            )
            exit()

        else:
            N_shm = int(sys.argv[1])

    case_file = '/path/to/caselist.txt'
    with open(case_file) as f:
        content = f.read()
    caselist = content.split()

    print('\n\n #### CONNECTOM site #### \n\n')

    print('========================================================================================')
    print('After bvalue mapping and resampling, volume difference')
    imgPath_py = '/path/to/dMRIharmonization/lib/tests/connectom_prisma/connectom/' \
        'XYZ/dwi_XYZ_connectom_st_b1200_resampled.nii.gz'
    imgPath_mat = '/path/to/dMRIharmonization/compare_matlab/connectom_prisma/connectom/' \
        'XYZ/rish/dwi_XYZ_connectom_st_b1200_bMapped_Resampled.nii.gz'
    dwi_difference(imgPath_mat, imgPath_py, caselist)
    print('========================================================================================')

    print('========================================================================================')
    print('After preprocessing, bval/bvec difference')
    # imgPath_mat and imgPath_py are nifti image paths from which bval/bvec
    # paths are understood
    bval_bvec_difference(imgPath_mat, imgPath_py, caselist)
    print('========================================================================================')

    print('========================================================================================')
    print('rish difference')
    imgPath_py = '/path/to/dMRIharmonization/lib/tests/connectom_prisma/connectom/' \
        'XYZ/harm/dwi_XYZ_connectom_st_b1200_resampled_LORDER.nii.gz'
    imgPath_mat = '/path/to/dMRIharmonization/compare_matlab/connectom_prisma/connectom/' \
        'XYZ/rish/rish/dwi_XYZ_connectom_st_b1200_bMapped_Resampled_LORDER.nii.gz'
    rish_diff(imgPath_mat, imgPath_py, caselist, N_shm)
    print('========================================================================================')

    print('========================================================================================')
    print('dti difference')
    imgPath_py = '/path/to/dMRIharmonization/lib/tests/connectom_prisma/connectom/' \
        'XYZ/dti/dwi_XYZ_connectom_st_b1200_resampled_FA.nii.gz'
    imgPath_mat = '/path/to/dMRIharmonization/compare_matlab/connectom_prisma/connectom/' \
        'XYZ/dti/dwi_XYZ_connectom_st_b1200_bMapped_Resampled_FA.nii.gz'
    dti_diff(imgPath_mat, imgPath_py, caselist)
    print('========================================================================================')


if __name__ == '__main__':
    main()
