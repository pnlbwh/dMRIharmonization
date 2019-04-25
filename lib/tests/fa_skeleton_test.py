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

from plumbum.cmd import antsApplyTransforms
from plumbum import FG
import multiprocessing
from numpy import mean
from test_util import *
import configparser, argparse

ROOTDIR= abspath(pjoin(LIBDIR, '..'))
mniTmp = pjoin(ROOTDIR, 'IITAtlas', 'IITmean_FA.nii.gz')

config = configparser.ConfigParser()
config.read(pjoin(LIBDIR,'config.ini'))
N_proc = int(config['DEFAULT']['N_proc'])
diffusionMeasures = [x for x in config['DEFAULT']['diffusionMeasures'].split(',')]

def read_caselist(file):

    with open(file) as f:

        imgs = []
        content= f.read()
        for line, row in enumerate(content.split()):
            temp= [element for element in row.split(',') if element] # handling w/space

            for img in temp:
                if not isfile(img):
                    raise FileNotFoundError(f'{img} does not exist: check line {line} in {file}')

            imgs.append(temp[0])

    return imgs

def antsReg(img, mask, mov, outPrefix, n_thread=1):

    if mask:
        check_call((' ').join(['antsRegistrationSyNQuick.sh',
                               '-d', '3',
                               '-f', img,
                               '-x', mask,
                               '-m', mov,
                               '-n', str(n_thread),
                               '-o', outPrefix]), shell= True)
    else:
        check_call((' ').join(['antsRegistrationSyNQuick.sh',
                               '-d', '3',
                               '-f', img,
                               '-m', mov,
                               '-n', str(n_thread),
                               '-o', outPrefix]), shell= True)


def register_subject(imgPath, warp2mni, trans2mni, templatePath, siteName):

    print(f'Warping {imgPath} diffusion measures to standard space')
    directory = dirname(imgPath)
    outPrefix = imgPath.split('.')[0] # should have _FA a the end
    prefix = psplit(outPrefix)[-1].replace('_FA', '')

    dmTmp = pjoin(templatePath, f'Mean_{siteName}_FA.nii.gz')
    maskTmp = pjoin(templatePath, f'{siteName}_Mask.nii.gz')
    warp2tmp = outPrefix + '1Warp.nii.gz'
    trans2tmp = outPrefix + '0GenericAffine.mat'
    # signal reconstruction might change with zero padding size, median filtering kernel size, and harmonized mask
    # so in case multiple debug is needed, redo the registration
    antsReg(dmTmp, maskTmp, imgPath, outPrefix)

    for dm in diffusionMeasures:
        output = pjoin(directory, prefix + f'_InMNI_{dm}.nii.gz')
        moving = pjoin(directory, prefix + f'_{dm}.nii.gz')
        # warp diffusion measure to template space first, then to MNI space
        antsApplyTransforms[
            '-d', '3',
            '-i', moving,
            '-o', output,
            '-r', mniTmp,
            '-t', warp2mni, trans2mni, warp2tmp, trans2tmp
        ] & FG

    return pjoin(directory, prefix + f'_InMNI_FA.nii.gz')

def sub2tmp2mni(templatePath, siteName, faImgs):

    # obtain the transform
    moving = pjoin(templatePath, f'Mean_{siteName}_FA.nii.gz')

    outPrefix= pjoin(templatePath, f'TemplateToMNI_{siteName}')
    warp2mni= outPrefix+'1Warp.nii.gz'
    trans2mni= outPrefix+'0GenericAffine.mat'
    # template is created once, it is expected that the user wants to keep the template same during debugging
    # so in case multiple debug is needed, pass the registration
    if not isfile(warp2mni):
        antsReg(mniTmp, None, moving, outPrefix, 8)


    pool= multiprocessing.Pool(N_proc)
    res=[]
    for imgPath in faImgs:
            res.append(pool.apply_async(func= register_subject,
                       args= (imgPath, warp2mni, trans2mni, templatePath, siteName, )))

    mniFAimgs= [r.get() for r in res]

    pool.close()
    pool.join()

    return mniFAimgs


def analyzeStat(faImgs):
    '''
    :param file: list of (FA or MD or GFA) that are already in MNI space
    :return: mean of the images
    '''

    skel= load(pjoin(ROOTDIR, 'IITAtlas', 'IITmean_FA_skeleton.nii.gz'))
    skel_mask= (skel.get_data()>0)*1.

    meanAttr=[]
    for faImg in faImgs:
        data= load(faImg).get_data()
        temp= data*skel_mask
        meanAttr.append(temp[temp>0].mean())

    return meanAttr


def main():

    parser = argparse.ArgumentParser(description='''Warps diffusion measures (FA, MD, GFA) to template space 
    and then to subject space. Finally, calculates mean FA over IITmean_FA_skeleton.nii.gz''')
    parser.add_argument('-i', '--input', type=str, required=True, help='input list of FA images')
    parser.add_argument('-s', '--site', type= str, required=True,
                        help='site name for locating template FA and mask in tempalte directory')
    parser.add_argument('-t', '--template', type=str, required=True,
                        help='template directory where Mean_{site}_FA.nii.gz and {site}_Mask.nii.gz is located')

    args = parser.parse_args()
    caselist=args.input
    siteName=args.site
    templatePath=args.template

    # read FA image list
    faImgs= read_caselist(caselist)

    # register and obtain *_InMNI_FA.nii.gz
    mniFAimgs= sub2tmp2mni(templatePath, siteName, faImgs)

    # pass *_InMNI_FA.nii.gz list to analyzeStat
    site_means= analyzeStat(mniFAimgs)
    print(f'{siteName} mean FA: ', mean(site_means))


if __name__ == '__main__':
    main()