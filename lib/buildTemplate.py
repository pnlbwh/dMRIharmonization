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

from plumbum.cmd import antsApplyTransforms
from plumbum import FG
from glob import glob

from scipy.ndimage import binary_opening, generate_binary_structure
from scipy.ndimage.filters import gaussian_filter
from util import *

eps= 2.2204e-16
SCRIPTDIR= os.path.dirname(__file__)
config = configparser.ConfigParser()
config.read(f'/tmp/harm_config_{os.getpid()}.ini')
N_shm = int(config['DEFAULT']['N_shm'])
N_proc = int(config['DEFAULT']['N_proc'])
diffusionMeasures= [x for x in config['DEFAULT']['diffusionMeasures'].split(',')]
travelHeads= int(config['DEFAULT']['travelHeads'])

def applyXform(inImg, refImg, warp, trans, outImg):

    antsApplyTransforms[
        '-d', '3',
        '-i', inImg,
        '-o', outImg,
        '-r', refImg,
        '-t', warp, '-t', trans
        ] & FG


def warp_bands(imgPath, maskPath, templatePath):

    prefix= os.path.basename(imgPath).split('.')[0]
    directory= os.path.dirname(imgPath)
    warp = glob(os.path.join(templatePath, prefix + f'_FA*[!Inverse]Warp.nii.gz'))
    trans = glob(os.path.join(templatePath, prefix + f'_FA*GenericAffine.mat'))

    # warping the mask
    applyXform(maskPath,
               os.path.join(templatePath, 'template0.nii.gz'),
               warp, trans,
               os.path.join(templatePath, os.path.abspath(maskPath).split('.')[0]+ 'Warped.nii.gz'))


    # warping the rish features
    for i in range(0, N_shm+1, 2):
        applyXform(os.path.join(directory, 'harm', f'{prefix}_L{i}.nii.gz'),
           os.path.join(templatePath, 'template0.nii.gz'),
           warp, trans,
           os.path.join(templatePath, f'{prefix}_WarpedL{i}.nii.gz'))


    # warping the diffusion measures
    for dm in diffusionMeasures:
        applyXform(os.path.join(directory, 'dti', f'{prefix}_{dm}.nii.gz'),
                   os.path.join(templatePath, 'template0.nii.gz'),
                   warp, trans,
                   os.path.join(templatePath, f'{prefix}_Warped{dm}.nii.gz'))


def createAntsCaselist(imgs, file):

    with open(file,'w') as f:
        for imgPath in imgs:
            prefix= os.path.basename(imgPath).split('.')[0]
            directory= os.path.dirname(imgPath)

            FA= os.path.join(directory,'dti', f'{prefix}_FA.nii.gz')
            L0= os.path.join(directory,'harm', f'{prefix}_L0.nii.gz')
            f.write(f'{FA},{L0}\n')


def antsMult(caselist, outPrefix):
    
    N_core=os.getenv('TEMPLATE_CONSTRUCT_CORES')
    check_call((' ').join([os.path.join(SCRIPTDIR, 'antsMultivariateTemplateConstruction2_fixed_random_seed.sh'),
                           '-d', '3',
                           '-g', '0.2',
                           '-k', '2',
                           '-t', "BSplineSyN[0.1,26,0]",
                           '-r', '1',
                           '-c', '2',
                           '-j', str(N_core) if N_core else str(N_proc),
                           '-f', '8x4x2x1',
                           '-o', outPrefix,
                           caselist]), shell= True)


def dti_stat(siteName, imgs, masks, templatePath, templateHdr):

    maskData = []
    for maskPath in masks:
        maskData.append(load_nifti(os.path.join(templatePath, os.path.abspath(maskPath).split('.')[0] + 'Warped.nii.gz'))[0])


    morphed_mask= binary_opening(np.mean(maskData, axis= 0)>0.5, structure= generate_binary_structure(3,1))*1
    morphed_mask_name= os.path.join(templatePath, f'{siteName}_Mask.nii.gz')
    templateAffine = templateHdr.get_best_affine()
    save_nifti(morphed_mask_name, morphed_mask.astype('uint8'), templateAffine, templateHdr)


    for dm in diffusionMeasures:
        imgData= []
        for imgPath in imgs:
            prefix = os.path.basename(imgPath).split('.')[0]
            imgData.append(load_nifti(os.path.join(templatePath, f'{prefix}_Warped{dm}.nii.gz'))[0])

        save_nifti(os.path.join(templatePath, f'Mean_{siteName}_{dm}.nii.gz'),
                                np.mean(imgData, axis= 0), templateAffine, templateHdr)

        save_nifti(os.path.join(templatePath, f'Std_{siteName}_{dm}.nii.gz'),
                                np.std(imgData, axis= 0), templateAffine, templateHdr)

    return morphed_mask_name


def rish_stat(siteName, imgs, templatePath, templateHdr):

    for i in range(0, N_shm+1, 2):
        imgData= []
        for imgPath in imgs:
            prefix = os.path.basename(imgPath).split('.')[0]
            imgData.append(load_nifti(os.path.join(templatePath, f'{prefix}_WarpedL{i}.nii.gz'))[0])

        templateAffine= templateHdr.get_best_affine()
        save_nifti(os.path.join(templatePath, f'Mean_{siteName}_L{i}.nii.gz'),
                            np.mean(imgData, axis= 0), templateAffine, templateHdr)

        save_nifti(os.path.join(templatePath, f'Std_{siteName}_L{i}.nii.gz'),
                            np.std(imgData, axis= 0), templateAffine, templateHdr)


def template_masking(refMaskPath, targetMaskPath, templatePath, siteName):

    ref= load(refMaskPath)
    target= load(targetMaskPath)

    templateMask= applymask(ref.get_data(), target.get_data())

    save_nifti(os.path.join(templatePath, 'templateMask.nii.gz'), templateMask.astype('uint8'), ref.affine, ref.header)

    for dm in diffusionMeasures:
        fileName= os.path.join(templatePath, f'Mean_{siteName}_{dm}.nii.gz')
        img= load(fileName)
        save_nifti(fileName, applymask(img.get_data(), templateMask), img.affine, img.header)

        fileName= os.path.join(templatePath, f'Std_{siteName}_{dm}.nii.gz')
        img= load(fileName)
        save_nifti(fileName, applymask(img.get_data(), templateMask), img.affine, img.header)


    return templateMask


def smooth(data):

    # Equivalence b/w a Gaussian and a Box filter
    # https://stackoverflow.com/questions/35340197/box-filter-size-in-relation-to-gaussian-filter-sigma
    # Matlab does [3x3x3] box smoothing by default
    sigma= np.repeat(np.sqrt(2/3),3)

    return gaussian_filter(data, sigma)


def stat_calc(ref, target, mask):

    ref= applymask(ref, mask)
    target= applymask(target, mask)

    delta= ref- target
    per_diff= 100*delta/(ref+eps)
    np.nan_to_num(per_diff).clip(max=100., min=-100., out= per_diff)
    per_diff_smooth= smooth(per_diff)
    scale= ref/(target+eps)

    return (delta, per_diff, per_diff_smooth, scale)


def difference_calc(refSite, targetSite, refImgs, targetImgs,
                    templatePath, templateHdr, mask, measures):

    '''
    if traveling heads:
        for each subject:
           delta= ref- target
        mean(delta)
        calc statistics
        save results
    else:
        delta= load(mean(ref))- load(mean(target))
        calc statistics
        save results
    '''

    templateAffine = templateHdr.get_best_affine()
    for dm in measures:
        delta=[]
        per_diff=[]
        per_diff_smooth= []
        scale= []
        if travelHeads:
            print('Using travelHeads for computing templates of',dm)
            for refImg, targetImg in zip(refImgs, targetImgs):
                prefix = os.path.basename(refImg).split('.')[0]
                ref= load_nifti(os.path.join(templatePath, f'{prefix}_Warped{dm}.nii.gz'))[0]

                prefix = os.path.basename(targetImg).split('.')[0]
                target= load_nifti(os.path.join(templatePath, f'{prefix}_Warped{dm}.nii.gz'))[0]

                temp= stat_calc(ref, target, mask)
                delta.append(temp[0])
                per_diff.append(temp[1])
                per_diff_smooth.append(temp[2])
                scale.append(temp[3])

        else:
            fileName= os.path.join(templatePath, f'Mean_{refSite}_{dm}.nii.gz')
            ref= load_nifti(fileName)[0]

            fileName= os.path.join(templatePath, f'Mean_{targetSite}_{dm}.nii.gz')
            target= load_nifti(fileName)[0]

            temp= stat_calc(ref, target, mask)
            delta.append(temp[0])
            per_diff.append(temp[1])
            per_diff_smooth.append(temp[2])
            scale.append(temp[3])

        save_nifti(os.path.join(templatePath, f'Delta_{dm}.nii.gz'),
                   np.mean(delta, axis= 0), templateAffine, templateHdr)

        save_nifti(os.path.join(templatePath, f'PercentageDiff_{dm}.nii.gz'),
                   np.mean(per_diff, axis= 0), templateAffine, templateHdr)

        save_nifti(os.path.join(templatePath, f'PercentageDiff_{dm}smooth.nii.gz'),
                   np.mean(per_diff_smooth, axis= 0), templateAffine, templateHdr)

        if 'L' in dm:
            save_nifti(os.path.join(templatePath, f'Scale_{dm}.nii.gz'),
                       np.sqrt(np.mean(scale, axis= 0)), templateAffine, templateHdr)



if __name__ == '__main__':
    pass
