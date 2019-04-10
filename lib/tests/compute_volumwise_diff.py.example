#!/usr/bin/env python

import warnings
with warnings.catch_warnings():
    warnings.filterwarnings("ignore", category=FutureWarning)
    import nibabel as nib

from conversion import read_bvals, read_bvecs
import numpy as np
import sys

np.set_printoptions(precision=5)
# BINS= [-np.inf, -10, -5, -1, 0, 1, 5, 10, np.inf]
BINS= [0, 1, 5, 10, np.inf]

def hist_calc(a, bins):

    a= np.array([item for sublist in a for item in sublist])

    hist_string= []
    N= len(bins)
    for i in range(N-1):
        hist_string.append(f'{bins[i]} <--> {bins[i+1]}')

    temp= np.histogram(a, bins=bins)[0]
    hist= temp.astype('float')/a.size

    try:
        hist[np.isnan]= 0
    except:
        pass

    print('%20s : %s' %('Bins','Density'))
    for i in range(N-1):
        print('%20s : %.6f' %(hist_string[i],hist[i]))

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
        bvals_py = np.array(read_bvals(imgPath_py.split('.')[0] + '.bval'))
        bvecs_py = np.array(read_bvecs(imgPath_py.split('.')[0] + '.bvec'))

        bvals_py, bvecs_py = stack_one_b0(bvals_py, bvecs= bvecs_py)


        # load matlab
        bvals_mat = np.array(read_bvals(imgPath_mat.split('.')[0] + '.bval'))
        bvecs_mat = np.array(read_bvecs(imgPath_mat.split('.')[0] + '.bvec'))

        # obtain case wise difference for bvals
        diff = abs(bvals_py.astype('float64') - bvals_mat.astype('float64'))
        bval_diff_list.append(diff.flatten())
        bval_total_diff.append(diff.sum())
        bval_mean_diff.append(diff.mean())

        # obtain case wise difference
        diff = abs(bvecs_py.astype('float64') - bvecs_mat.astype('float64'))
        bvec_diff_list.append(diff.flatten())
        bvec_total_diff.append(diff.sum())
        bvec_mean_diff.append(diff.mean())

    print(f'Mean of total difference over all the bvals: {np.mean(bval_total_diff)}')
    print(f'Mean of mean difference over all the bvals: {np.mean(bval_mean_diff)}')
    hist_calc(bval_diff_list, bins= BINS)
    print('\n\n')

    print(f'Mean of total difference over all the bvecs: {np.mean(bvec_total_diff)}')
    print(f'Mean of mean difference over all the bvecs: {np.mean(bvec_mean_diff)}')
    hist_calc(bvec_diff_list, bins= BINS)
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
        img_py = nib.load(imgPath_py).get_data()

        # keep only the first b0
        bvals = np.array(read_bvals(imgPath_py.split('.')[0] + '.bval'))
        img_py= stack_one_b0(bvals, dwi= img_py)

        pyX, pyY, pyZ, _ = img_py.shape

        # load matlab
        img_mat = nib.load(imgPath_mat).get_data()[:pyX, :pyY, :pyZ,: ]

        # obtain case wise difference
        diff = abs(img_py.astype('float64') - img_mat.astype('float64'))
        diff_list.append(diff.flatten())
        total_diff.append(diff.sum())
        mean_diff.append(diff.mean())

        # save the difference mask
        affine = nib.load(imgPath_py).affine
        nib.save(nib.Nifti1Image(diff, affine), imgPath_mat.split('.')[0]+'_diff_mask.nii.gz')



    print(f'Mean of total difference over all the voxels: {np.mean(total_diff)}')
    print(f'Mean of mean difference over all the voxels: {np.mean(mean_diff)}')
    hist_calc(diff_list, bins= BINS)


def stack_one_b0(bvals, bvecs= None, dwi=None):

    bvals= np.array(bvals)
    where_b0 = np.where(bvals == 0)[0]
    where_dwi = np.where(bvals != 0)[0]

    reduced_ind= [where_b0[0]]
    for i in where_dwi:
        reduced_ind.append(i)


    if dwi is not None:
        reduced_volume= np.take(dwi, indices=reduced_ind, axis=3)
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
        img_py = nib.load(imgPath_py).get_data()
        pyX, pyY, pyZ = img_py.shape

        # load matlab
        img_mat = nib.load(imgPath_mat).get_data()[:pyX, :pyY, :pyZ]

        # obtain case wise difference
        diff = abs(img_py.astype('float64') - img_mat.astype('float64'))
        diff_list.append(diff.flatten())
        total_diff.append(diff.sum())
        mean_diff.append(diff.mean())

    print(f'Mean of total difference over all the voxels: {np.mean(total_diff)}')
    print(f'Mean of mean difference over all the voxels: {np.mean(mean_diff)}')
    hist_calc(diff_list, bins= BINS)
    print('\n\n')

# rish difference
def rish_diff(imgPath_given_mat, imgPath_given_py, caselist, N_shm= 6):

    for i in range(0, N_shm + 1, 2):

        imgPath_mat = imgPath_given_mat.replace('ORDER', str(i))
        imgPath_py = imgPath_given_py.replace('ORDER', str(i))

        diff_list = []
        total_diff= []
        mean_diff= []

        for case in caselist:
            imgPath_mat = imgPath_mat.replace('XYZ', case)
            imgPath_py = imgPath_py.replace('XYZ', case)

            # load python
            img_py = nib.load(imgPath_py).get_data()
            pyX, pyY, pyZ = img_py.shape

            # load matlab
            img_mat = nib.load(imgPath_mat).get_data()[:pyX, :pyY, :pyZ]

            # obtain case wise difference
            diff = abs(img_py.astype('float64') - img_mat.astype('float64'))
            diff_list.append(diff.flatten())
            total_diff.append(diff.sum())
            mean_diff.append(diff.mean())

        print(f'Feature: L{i}: ')
        print(f'Mean of total difference over all the voxels: {np.mean(total_diff)}')
        print(f'Mean of mean difference over all the voxels: {np.mean(mean_diff)}')
        hist_calc(diff_list, bins= BINS)
        print('\n\n')


# scale difference
def scale_diff(imgPath_given_mat, imgPath_given_py, N_shm= 6):

    for i in range(0, N_shm + 1, 2):

        imgPath_mat = imgPath_given_mat.replace('ORDER', str(i))
        imgPath_py = imgPath_given_py.replace('ORDER', str(i))

        diff_list = []
        total_diff= []
        mean_diff= []

        # load python
        img_py = nib.load(imgPath_py).get_data()
        pyX, pyY, pyZ = img_py.shape

        # load matlab
        img_mat = nib.load(imgPath_mat).get_data()[:pyX, :pyY, :pyZ]

        # obtain case wise difference
        diff = abs(img_py.astype('float64') - img_mat.astype('float64'))
        diff_list.append(diff.flatten())
        total_diff.append(diff.sum())
        mean_diff.append(diff.mean())

        print(f'Feature: L{i}: ')
        print(f'Mean of total difference over all the voxels: {np.mean(total_diff)}')
        print(f'Mean of mean difference over all the voxels: {np.mean(mean_diff)}')
        hist_calc(diff_list, bins= BINS)
        print('\n\n')


def main():
    if sys.argv[1]=='-h' or sys.argv[1]=='--help':
        print(
'''
This module is the gateway for testing equivalence between MATALAB and PYTHON harmonization results. 
Edit 'compute_volumwise_diff.py' the module for computing volumewise difference over all cases in a site.
Specify a sample path in 'imgPath_py' and 'imgPath_mat' using *XYZ* for case, *ORDER* for shm order in sample path.
The program will replace *XYZ* with caseid obtained from provided caselist. It will also substitute *ORDER* with
proper spherical harmonic order upto N_shm given to
rish_diff(imgPath_given_mat, imgPath_given_py, caselist, N_shm= 6) and
scale_diff(imgPath_given_mat, imgPath_given_py, N_shm= 6)
'''
            )
        exit()

    case_file = '.../Harmonization-Python/compare_matlab/BSNIP_Baltimore/CIDAR-post/' \
                'caselist.txt'
    with open(case_file) as f:
        content= f.read()
    caselist= content.split()



    print('#### BSNIP_Baltimore, reference site: CIDAR-post #### \n\n')

    print('========================================================================================')
    imgPath_py='.../Harmonization-Python/BSNIP_Baltimore/CIDAR-post/' \
                'caseXYZ/XYZ-dwi-Ed-centered_bMapped.nii.gz'
    imgPath_mat='.../Harmonization-Python/compare_matlab/BSNIP_Baltimore/CIDAR-post/' \
                'caseXYZ/XYZ-dwi-Ed-centered_bMapped.nii.gz'
    print('After bvalue mapping, volume difference')
    dwi_difference(imgPath_mat, imgPath_py, caselist)

    print('========================================================================================')
    imgPath_py='.../Harmonization-Python/BSNIP_Baltimore/CIDAR-post/' \
                'caseXYZ/XYZ-dwi-Ed-centered_resampled.nii.gz'
    imgPath_mat='.../Harmonization-Python/compare_matlab/BSNIP_Baltimore/CIDAR-post/' \
                'caseXYZ/XYZ-dwi-Ed-centered_Resampled.nii.gz'
    print('After resampling, volume difference')
    dwi_difference(imgPath_mat, imgPath_py, caselist)

    print('========================================================================================')
    imgPath_py='.../Harmonization-Python/BSNIP_Baltimore/CIDAR-post/' \
                'caseXYZ/XYZ-dwi-Ed-centered_resampled.nii.gz'
    imgPath_mat='.../Harmonization-Python/compare_matlab/BSNIP_Baltimore/CIDAR-post/' \
                'caseXYZ/rish/XYZ-dwi-Ed-centered_bMapped_Resampled.nii.gz'
    print('After preprocessing, volume difference')
    dwi_difference(imgPath_mat, imgPath_py, caselist)

    print('========================================================================================')
    print('After preprocessing, bval/bvec difference')
    bval_bvec_difference(imgPath_mat, imgPath_py, caselist)

    print('========================================================================================')
    print('After preprocessing, mask difference')
    imgPath_py='.../Harmonization-Python/BSNIP_Baltimore/CIDAR-post/' \
                'caseXYZ/Tensor_mask-XYZ-dwi-filt-Ed_AvGradient-cleaned_resampled.nii.gz'
    imgPath_mat='.../Harmonization-Python/compare_matlab/BSNIP_Baltimore/CIDAR-post/' \
                'caseXYZ/rish/Tensor_mask-XYZ-dwi-filt-Ed_AvGradient-cleaned_Resampled.nii.gz'
    dti_diff(imgPath_mat, imgPath_py, caselist)

    print('========================================================================================')
    print('After preprocessing, Rish difference')
    imgPath_py='.../Harmonization-Python/BSNIP_Baltimore/CIDAR-post/' \
                'caseXYZ/harm/XYZ-dwi-Ed-centered_resampled_LORDER.nii.gz'
    imgPath_mat='.../Harmonization-Python/compare_matlab/BSNIP_Baltimore/CIDAR-post/' \
                'caseXYZ/rish/rish/XYZ-dwi-Ed-centered_bMapped_Resampled_LORDER.nii.gz'
    rish_diff(imgPath_mat, imgPath_py, caselist)

    print('========================================================================================')
    print('After preprocessing, DTI difference')
    imgPath_py='.../Harmonization-Python/BSNIP_Baltimore/CIDAR-post/' \
                'caseXYZ/dti/XYZ-dwi-Ed-centered_resampled_FA.nii.gz'
    imgPath_mat='.../Harmonization-Python/compare_matlab/BSNIP_Baltimore/CIDAR-post/' \
                'caseXYZ/dti/XYZ-dwi-Ed-centered_bMapped_Resampled_FA.nii.gz'
    dti_diff(imgPath_mat, imgPath_py, caselist)

    # ========================================================================================================================


    print('#### Template for BSNIP and CIDAR #### \n')

    print('========================================================================================')
    print('After template creation, rish scale map difference')
    imgPath_py='.../Harmonization-Python/BSNIP_Baltimore/template/Scale_LORDER.nii.gz'
    imgPath_mat='.../Harmonization-Python/compare_matlab/BSNIP_Baltimore/template/Scale_LORDER.nii.gz'
    scale_diff(imgPath_mat, imgPath_py)


    # ========================================================================================================================

    case_file = '.../Harmonization-Python/compare_matlab/BSNIP_Baltimore/BSNIP_Balt_trainingHC/' \
                'caselist.txt'
    with open(case_file) as f:
        content= f.read()
    caselist= content.split()

    print('#### BSNIP_Baltimore, target site: BSNIP #### \n')

    print('========================================================================================')
    imgPath_py='.../Harmonization-Python/BSNIP_Baltimore/BSNIP_Balt_trainingHC/' \
                'XYZ/XYZ_dwi_xc_Ed_resampled.nii.gz'
    imgPath_mat='.../Harmonization-Python/compare_matlab/BSNIP_Baltimore/BSNIP_Balt_trainingHC/' \
                'XYZ/rish/XYZ_dwi_xc_Ed_Resampled.nii.gz'
    print('After preprocessing, before harmonization, volume difference')
    dwi_difference(imgPath_mat, imgPath_py, caselist)

    print('========================================================================================')
    print('After preprocessing, before harmonization, bval/bvec difference')
    bval_bvec_difference(imgPath_mat, imgPath_py, caselist)

    print('========================================================================================')
    print('After preprocessing, before harmonization mask difference')
    imgPath_py='.../Harmonization-Python/BSNIP_Baltimore/BSNIP_Balt_trainingHC/' \
                'XYZ/XYZ_dwi_xc_Ed_OTSUtensormask_cleaned_resampled.nii.gz'
    imgPath_mat='.../Harmonization-Python/compare_matlab/BSNIP_Baltimore/BSNIP_Balt_trainingHC/' \
                'XYZ/rish/XYZ_dwi_xc_Ed_OTSUtensormask_cleaned_Resampled.nii.gz'
    dti_diff(imgPath_mat, imgPath_py, caselist)

    print('========================================================================================')
    print('After preprocessing, before harmonization, rish difference')
    imgPath_py='.../Harmonization-Python/BSNIP_Baltimore/BSNIP_Balt_trainingHC/' \
                'XYZ/harm/XYZ_dwi_xc_Ed_resampled_LORDER.nii.gz'
    imgPath_mat='.../Harmonization-Python/compare_matlab/BSNIP_Baltimore/BSNIP_Balt_trainingHC/' \
                'XYZ/rish/rish/XYZ_dwi_xc_Ed_Resampled_LORDER.nii.gz'
    rish_diff(imgPath_mat, imgPath_py, caselist)

    print('========================================================================================')
    print('After preprocessing, before harmonization, DTI difference')
    imgPath_py='.../Harmonization-Python/BSNIP_Baltimore/BSNIP_Balt_trainingHC/' \
                'XYZ/dti/XYZ_dwi_xc_Ed_resampled_FA.nii.gz'
    imgPath_mat='.../Harmonization-Python/compare_matlab/BSNIP_Baltimore/BSNIP_Balt_trainingHC/' \
                'XYZ/dti/XYZ_dwi_xc_Ed_Resampled_FA.nii.gz'
    dti_diff(imgPath_mat, imgPath_py, caselist)


    print('========================================================================================')
    imgPath_py='.../Harmonization-Python/BSNIP_Baltimore/BSNIP_Balt_trainingHC/' \
                'XYZ/harmonized_XYZ_dwi_xc_Ed_resampled.nii.gz'
    imgPath_mat='.../Harmonization-Python/compare_matlab/BSNIP_Baltimore/BSNIP_Balt_trainingHC/' \
                'XYZ/rish/harmonized_XYZ_dwi_xc_Ed_Resampled.nii.gz'
    print('After harmonization, volume difference \n')
    dwi_difference(imgPath_mat, imgPath_py, caselist,
                   '.../Harmonization-Python/compare_matlab/BSNIP_Baltimore/CIDAR-post/diff_mask_reconst.nii.gz')
    print('After harmonization, bval/bvec difference \n')
    bval_bvec_difference(imgPath_mat, imgPath_py, caselist)

    print('========================================================================================')
    print('After harmonization, mask difference: ')
    imgPath_py='.../Harmonization-Python/BSNIP_Baltimore/BSNIP_Balt_trainingHC/' \
                'XYZ/harmonized_XYZ_dwi_xc_Ed_resampled_mask.nii.gz'
    imgPath_mat='.../Harmonization-Python/compare_matlab/BSNIP_Baltimore/BSNIP_Balt_trainingHC/' \
                'XYZ/rish/harmonized_XYZ_dwi_xc_Ed_OTSUtensormask_cleaned_Resampled.nii.gz'
    dti_diff(imgPath_mat, imgPath_py, caselist)



if __name__=='__main__':
    main()
