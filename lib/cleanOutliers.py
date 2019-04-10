from subprocess import check_call
import os, configparser, shutil
from skimage.measure import label, regionprops
from scipy.ndimage import binary_erosion, binary_dilation, \
    generate_binary_structure, iterate_structure
import numpy as np

import warnings
with warnings.catch_warnings():
    warnings.filterwarnings("ignore", category=FutureWarning)
    from dipy.io.image import load_nifti, save_nifti
    from dipy.segment.mask import applymask

from buildTemplate import applyXform
from rish import rish
from local_med_filter import local_med_filter
from preprocess import dti_harm, preprocessing

eps= 2.2204e-16
SCRIPTDIR= os.path.dirname(__file__)
config = configparser.ConfigParser()
config.read(os.path.join(SCRIPTDIR,'config.ini'))
N_shm = int(config['DEFAULT']['N_shm'])
N_proc = int(config['DEFAULT']['N_proc'])
debug = int(config['DEFAULT']['debug'])

def antsReg(img, mask, mov, outPrefix):

    if mask:
        check_call((' ').join(['antsRegistrationSyNQuick.sh',
                               '-d', '3',
                               '-f', img,
                               '-x', mask,
                               '-m', mov,
                               '-o', outPrefix]), shell= True)
    else:
        check_call((' ').join(['antsRegistrationSyNQuick.sh',
                               '-d', '3',
                               '-f', img,
                               '-m', mov,
                               '-o', outPrefix]), shell= True)

def antsApply(templatePath, directory, prefix):

    for i in range(0, N_shm+1, 2):
        moving= os.path.join(templatePath, f'Scale_L{i}.nii.gz')
        fixed= os.path.join(directory, f'{prefix}_L{i}.nii.gz')
        output= os.path.join(directory, f'Scale_L{i}_{prefix}.nii.gz')
        warp= os.path.join(directory, f'ToSubjectSpace_{prefix}1Warp.nii.gz')
        trans= os.path.join(directory, f'ToSubjectSpace_{prefix}0GenericAffine.mat')

        applyXform(moving, fixed, warp, trans, output)


def custom_spherical_structure(D):
    '''Given diameter D, returns a spherical structuring element'''

    sw=(D-1)/2
    ses2= int(np.ceil(D/2))
    [y,x,z]= np.meshgrid(np.arange(-sw,sw+1),np.arange(-sw,sw+1),np.arange(-sw,sw+1))
    m= np.sqrt(x**2 + y**2 + z**2)
    b= (m <= m[ses2,ses2,D-1])

    return b


def findLargestConnectMask(img, mask):

    mask_scale = label(img > 0.00001, connectivity=1)
    maxArea = 0
    for region in regionprops(mask_scale):
        if region.area > maxArea:
            maxLabel = region.label
            maxArea = region.area

    largeConnectMask = (mask_scale == maxLabel)
    mask *= largeConnectMask

    return mask


def ring_masking(directory, prefix, maskPath, shm_coeff, b0, qb_model):

    B = qb_model.B
    bvals= qb_model.gtab.bvals

    mapped_cs= []
    shs_same_level= [[0, 1], [1, 6], [6, 15], [15, 28], [28, 45]]

    mask, _ = load_nifti(maskPath)
    for i in range(0, N_shm+1, 2):

        # load data and mask
        fileName= os.path.join(directory, 'harm', f'Scale_L{i}_{prefix}.nii.gz')
        img, affine = load_nifti(fileName)

        if i==0: # compute the skullRingMask from 0th shm

            mask= findLargestConnectMask(img, mask)

            n_zero= 5
            se= custom_spherical_structure(n_zero)
            paddedMask= np.pad(mask, n_zero, 'constant', constant_values= 0.)

            dilM = binary_dilation(paddedMask, se)*1
            eroM = binary_erosion(paddedMask, se)*1
            skullRingMask = dilM - eroM

        paddedImg = np.pad(img, n_zero, 'constant', constant_values=0.)
        skullRing= paddedImg*skullRingMask
        thresh= np.percentile(skullRing[skullRing>0], 95)
        outLier= (skullRing>thresh)*1
        tmp= local_med_filter(paddedImg, outLier)
        denoisedImg = tmp[n_zero:-n_zero, n_zero:-n_zero, n_zero:-n_zero]


        # for i=0, create new mask over denoisedImg, save it as harmonized_{prefix}_mask
        if i==0:
            mask_final= findLargestConnectMask(denoisedImg, mask)
            # mask_final= binary_dilation(mask_final, generate_binary_structure(3,1))*1
            harmMask = os.path.join(directory, f'harmonized_{prefix}_mask.nii.gz')
            save_nifti(harmMask, mask_final, affine=affine)


        ind= int(i/2)
        denoisedImg= applymask(denoisedImg, mask_final)
        for level in range(shs_same_level[ind][0], shs_same_level[ind][1]):
            mapped_cs.append(denoisedImg * shm_coeff[ :,:,:,level])

    S_hat= np.dot(np.moveaxis(mapped_cs, 0, -1), B.T)
    np.nan_to_num(S_hat).clip(min= 0., max= 1., out= S_hat)


    # affine= templateAffine for all Scale_L{i}
    mappedFile= os.path.join(directory, f'{prefix}_mapped_cs.nii.gz')
    save_nifti(mappedFile, S_hat, affine= affine)

    # un-normalize harmonized data
    S_hat_dwi= applymask(S_hat, b0) # overriding applymask function with a nonbinary mask b0

    # place b0s in proper indices
    S_hat_final= stack_b0(bvals, S_hat_dwi, b0)
    S_hat_final= applymask(S_hat_final, mask_final)


    # save harmonized data
    harmImg= os.path.join(directory, f'harmonized_{prefix}.nii.gz')
    save_nifti(harmImg, S_hat_final, affine= affine)


    return (harmImg, harmMask)


def reconst(imgPath, maskPath, moving, templatePath, preFlag):

    if preFlag:
        imgPath, maskPath = preprocessing(imgPath, maskPath)

    b0, shm_coeff, qb_model = dti_harm(imgPath, maskPath)

    directory = os.path.dirname(imgPath)
    inPrefix = imgPath.split('.')[0]
    prefix = os.path.split(inPrefix)[-1]

    print(f'Registering template FA to {imgPath} space ...')
    outPrefix = os.path.join(directory, 'harm', 'ToSubjectSpace_' + prefix)
    fixed = os.path.join(directory, 'dti', f'{prefix}_FA.nii.gz')
    antsReg(fixed, maskPath, moving, outPrefix)
    antsApply(templatePath, os.path.join(directory, 'harm'), prefix)

    print(f'Reconstructing signal from {imgPath} rish features ...')
    harmImg, harmMask = ring_masking(directory, prefix, maskPath, shm_coeff, b0, qb_model)
    shutil.copyfile(inPrefix + '.bvec', harmImg.split('.')[0] + '.bvec')
    shutil.copyfile(inPrefix + '.bval', harmImg.split('.')[0] + '.bval')

    if debug:
        dti_harm(harmImg, harmMask)

    return (imgPath, maskPath, harmImg, harmMask)

def stack_b0(bvals, dwi, b0):

    ind= np.where(bvals==0)[0]

    S_hat_final= []
    j= 0
    for i in range(len(bvals)):
        if i in ind:
            S_hat_final.append(b0)
        else:
            S_hat_final.append(dwi[:,:,:,j])
            j+=1

    return np.moveaxis(np.array(S_hat_final), 0, -1)

if __name__ =='__main__':
    pass