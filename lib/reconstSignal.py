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

from skimage.measure import label, regionprops
from scipy.ndimage import binary_erosion, binary_dilation, generate_binary_structure, iterate_structure

from util import *
from buildTemplate import applyXform
from local_med_filter import local_med_filter
from preprocess import dti_harm, preprocessing
from rish import rish

eps= 2.2204e-16
SCRIPTDIR= os.path.dirname(__file__)
config = configparser.ConfigParser()
config.read(f'/tmp/harm_config_{os.getpid()}.ini')
N_shm = int(config['DEFAULT']['N_shm'])
N_proc = int(config['DEFAULT']['N_proc'])
debug = int(config['DEFAULT']['debug'])
n_zero = int(config['DEFAULT']['N_zero'])

def antsReg(img, mask, mov, outPrefix):

    if mask:
        p= Popen((' ').join(['antsRegistrationSyNQuick.sh',
                               '-d', '3',
                               '-f', img,
                               '-x', mask,
                               '-m', mov,
                               '-o', outPrefix,
                               '-e', '123456']), shell= True)
        p.wait()
    else:
        p= Popen((' ').join(['antsRegistrationSyNQuick.sh',
                               '-d', '3',
                               '-f', img,
                               '-m', mov,
                               '-o', outPrefix,
                               '-e', '123456']), shell= True)
        p.wait()

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


def ring_masking(directory, prefix, maskPath, shm_coeff, b0, qb_model, hdr):

    B = qb_model.B
    affine= hdr.get_best_affine()

    mapped_cs= []
    shs_same_level= [[0, 1], [1, 6], [6, 15], [15, 28], [28, 45]]

    mask= load(maskPath).get_data()
    for i in range(0, N_shm+1, 2):

        # load data and mask
        fileName= os.path.join(directory, 'harm', f'Scale_L{i}_{prefix}.nii.gz')
        img= load(fileName).get_data()

        if i==0: # compute the skullRingMask from 0th shm

            mask= findLargestConnectMask(img, mask)

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
            save_nifti(harmMask, mask_final.astype('uint8'), affine, hdr)


        ind= int(i/2)
        denoisedImg= applymask(denoisedImg, mask_final)
        for level in range(shs_same_level[ind][0], shs_same_level[ind][1]):
            mapped_cs.append(denoisedImg * shm_coeff[ :,:,:,level])

    S_hat= np.dot(np.moveaxis(mapped_cs, 0, -1), B.T)
    # keep only upper half of the reconstructed signal
    S_hat= S_hat[..., :int(S_hat.shape[3]/2)]
    np.nan_to_num(S_hat).clip(min= 0., max= 1., out= S_hat)

    # affine= templateAffine for all Scale_L{i}
    mappedFile= os.path.join(directory, f'{prefix}_mapped_cs.nii.gz')
    save_nifti(mappedFile, S_hat, affine, hdr)

    # un-normalize harmonized data
    S_hat_dwi= applymask(S_hat, b0) # overriding applymask function with a nonbinary mask b0

    # place b0s in proper indices
    S_hat_final= stack_b0(qb_model.gtab.b0s_mask, S_hat_dwi, b0)
    S_hat_final= applymask(S_hat_final, mask_final)


    # save harmonized data
    harmImg= os.path.join(directory, f'harmonized_{prefix}.nii.gz')
    save_nifti(harmImg, S_hat_final, affine, hdr)


    return (harmImg, harmMask)


def reconst(imgPath, maskPath, moving, templatePath, preFlag):

    if preFlag:
        imgPath, maskPath = preprocessing(imgPath, maskPath)

    img = load(imgPath)

    directory = os.path.dirname(imgPath)
    inPrefix = imgPath.split('.nii')[0]
    prefix = os.path.split(inPrefix)[-1] 
    outPrefix = os.path.join(directory, 'harm', prefix) 
    b0, shm_coeff, qb_model = rish(imgPath, maskPath, inPrefix, outPrefix, N_shm)


    print(f'Registering template FA to {imgPath} space ...')
    outPrefix = os.path.join(directory, 'harm', 'ToSubjectSpace_' + prefix)
    fixed = os.path.join(directory, 'dti', f'{prefix}_FA.nii.gz')
    antsReg(fixed, maskPath, moving, outPrefix)
    antsApply(templatePath, os.path.join(directory, 'harm'), prefix)

    print(f'Reconstructing signal from {imgPath} rish features ...')
    harmImg, harmMask = ring_masking(directory, prefix, maskPath, shm_coeff, b0, qb_model, img.header)
    shutil.copyfile(inPrefix + '.bvec', harmImg.split('.nii')[0] + '.bvec')
    shutil.copyfile(inPrefix + '.bval', harmImg.split('.nii')[0] + '.bval')

    return (imgPath, maskPath, harmImg, harmMask)


def stack_b0(b0s_mask, dwi, b0):

    N= int(len(b0s_mask)/2)
    ind= np.where(b0s_mask[ :N])[0]

    S_hat_final= []
    j= 0
    for i in range(N):
        if i in ind:
            S_hat_final.append(b0)
        else:
            S_hat_final.append(dwi[:,:,:,j])
            j+=1

    return np.moveaxis(np.array(S_hat_final), 0, -1)
