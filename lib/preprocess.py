import os, shutil, configparser, multiprocessing
import numpy as np
from conversion import nifti_write


import warnings
with warnings.catch_warnings():
    warnings.filterwarnings("ignore", category=FutureWarning)
    from dipy.io.image import load_nifti, save_nifti
    from dipy.segment.mask import applymask
    from dipy.io import read_bvals_bvecs
    from nibabel import load

from .denoising import denoising
from .bvalMap import remapBval
from .resampling import resampling

from .dti import dti
from .rish import rish

SCRIPTDIR= os.path.dirname(__file__)
config = configparser.ConfigParser()
config.read(os.path.join(SCRIPTDIR,'config.ini'))

N_shm = int(config['DEFAULT']['N_shm'])
N_proc = int(config['DEFAULT']['N_proc'])
denoise= int(config['DEFAULT']['denoise'])
bvalMap= float(config['DEFAULT']['bvalMap'])
resample= config['DEFAULT']['resample']
if resample=='0':
    resample = 0

def write_bvals(bval_file, bvals):
    with open(bval_file, 'w') as f:
        f.write(('\n').join(str(b) for b in bvals))

def read_caselist(file):

    with open(file) as f:

        imgs = []
        masks = []
        content= f.read()
        for line, row in enumerate(content.split()):
            temp= [element for element in row.split(',') if element] # handling w/space
            imgs.append(temp[0])
            masks.append(temp[1])


        return (imgs, masks)


def dti_harm(imgPath, maskPath):

    directory = os.path.dirname(imgPath)
    inPrefix = imgPath.split('.')[0]
    prefix = os.path.split(inPrefix)[-1]

    outPrefix = os.path.join(directory, 'dti', prefix)

    # if the dti output exists with the same prefix, don't dtifit again
    if not os.path.exists(outPrefix+'_FA.nii.gz'):
        dti(imgPath, maskPath, inPrefix, outPrefix)

    outPrefix = os.path.join(directory, 'harm', prefix)
    b0, shm_coeff, qb_model= rish(imgPath, maskPath, inPrefix, outPrefix, N_shm)

    return (b0, shm_coeff, qb_model)


# def pre_dti_harm(imgPath, maskPath):
def pre_dti_harm(itr):
    imgPath, maskPath = preprocessing(itr[0], itr[1])
    dti_harm(imgPath, maskPath)
    return (imgPath, maskPath)

# convert NRRD to NIFTI on the fly
def nrrd2nifti(imgPath):

    if imgPath.endswith('.nrrd') or imgPath.endswith('.nhdr'):
        niftiImgPrefix= imgPath.split('.')[0]
        nifti_write(imgPath, niftiImgPrefix)

        return niftiImgPrefix+'.nii.gz'


def preprocessing(imgPath, maskPath):

    # load signal attributes for pre-processing ----------------------------------------------------------------
    imgPath= nrrd2nifti(imgPath)
    lowRes = load(imgPath)
    lowResImg = lowRes.get_data()
    lowResImgHdr = lowRes.header

    maskPath= nrrd2nifti(maskPath)
    lowRes = load(maskPath)
    lowResMask = lowRes.get_data()
    lowResMaskHdr = lowRes.header

    lowResImg = applymask(lowResImg, lowResMask)

    inPrefix = imgPath.split('.')[0]

    bvals, _ = read_bvals_bvecs(inPrefix + '.bval', None)

    # pre-processing -------------------------------------------------------------------------------------------
    suffix= None
    # modifies data only
    if denoise:
        print('Denoising ', imgPath)
        lowResImg, _ = denoising(lowResImg, lowResMask)
        suffix = '_denoised'
        # save_nifti(imgPath.split('.')[0]+'_denoised.nii.gz', lowResImg, lowResImgHdr.affine)

    # modifies data, and bvals
    if bvalMap:
        print('B value mapping ', imgPath)
        lowResImg, bvals = remapBval(lowResImg, bvals, bvalMap)
        suffix = '_bmapped'
        # save_nifti(imgPath.split('.')[0]+'_bmapped.nii.gz', lowResImg, lowResImgHdr.affine)

    # modifies data, mask, and headers
    if resample:
        print('Resampling ', imgPath)
        sp_high = np.array([float(i) for i in resample.split('x')])
        if (abs(sp_high-lowResImgHdr['pixdim'][1:4])>10e-3).any():
            imgPath, maskPath = \
                resampling(imgPath, maskPath, lowResImg, lowResImgHdr, lowResMask, lowResMaskHdr, sp_high, bvals)
            suffix = '_resampled'


    # save pre-processed data; resampled data is saved inside resampling() -------------------------------------
    if (denoise or bvalMap) and suffix!= '_resampled':
        imgPath = inPrefix + suffix + '.nii.gz'
        save_nifti(imgPath, lowResImg, lowResImgHdr.get_qform())

    if suffix:
        shutil.copyfile(inPrefix + '.bvec', inPrefix + suffix + '.bvec')
    if bvalMap:
        write_bvals(inPrefix + suffix + '.bval', bvals)
    elif denoise or suffix== '_resampled':
        shutil.copyfile(inPrefix + '.bval', inPrefix + suffix + '.bval')


    return (imgPath, maskPath)


def common_processing(caselist):
    imgs, masks = read_caselist(caselist)
    f = open(caselist + '.modified', 'w')

    pool = multiprocessing.Pool(N_proc)  # Use all available cores, otherwise specify the number you want as an argument

    res = pool.map_async(pre_dti_harm, np.hstack((np.reshape(imgs, (len(imgs), 1)), np.reshape(masks, (len(masks), 1)))))
    attributes = res.get()
    for i in range(len(imgs)):
        imgs[i] = attributes[i][0]
        masks[i] = attributes[i][1]
        f.write(f'{imgs[i]},{masks[i]}\n')

    pool.close()
    pool.join()

    f.close()

    return (imgs, masks)