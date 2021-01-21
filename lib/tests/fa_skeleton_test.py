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
import numpy as np
from test_util import *
import argparse
from conversion import read_imgs, read_imgs_masks
from harm_plot import harm_plot, generate_csv

ROOTDIR= abspath(pjoin(LIBDIR, '..'))
mniTmp = pjoin(ROOTDIR, 'IITAtlas', 'IITmean_FA.nii.gz')

diffusionMeasures = ['MD', 'FA', 'GFA']

def antsReg(img, mask, mov, outPrefix, n_thread=1):

    if mask:
        p= Popen((' ').join(['antsRegistrationSyNQuick.sh',
                               '-d', '3',
                               '-f', img,
                               '-x', mask,
                               '-m', mov,
                               '-n', str(n_thread),
                               '-o', outPrefix]), shell= True)
        p.wait()
    else:
        p= Popen((' ').join(['antsRegistrationSyNQuick.sh',
                               '-d', '3',
                               '-f', img,
                               '-m', mov,
                               '-n', str(n_thread),
                               '-o', outPrefix]), shell= True)
        p.wait()


def register_subject(imgPath, warp2mni, trans2mni, templatePath, siteName):

    print(f'Warping {imgPath} diffusion measures to standard space')
    directory = dirname(imgPath)
    outPrefix = pjoin(templatePath, imgPath.split('.nii')[0]) # should have _FA at the end
    prefix = psplit(outPrefix)[-1].replace('_FA', '')
    
    dmTmp = pjoin(templatePath, f'Mean_{siteName}_FA.nii.gz')
    maskTmp = pjoin(templatePath, f'{siteName}_Mask.nii.gz')
    warp2tmp = outPrefix + '1Warp.nii.gz'
    trans2tmp = outPrefix + '0GenericAffine.mat'
    # signal reconstruction might change with zero padding size, median filtering kernel size, and harmonized mask
    # so in case multiple debug is needed, redo the registration
    antsReg(dmTmp, maskTmp, imgPath, outPrefix)

    for dm in diffusionMeasures:
        output = pjoin(templatePath, prefix + f'_InMNI_{dm}.nii.gz')
        moving = pjoin(directory, prefix + f'_{dm}.nii.gz')
        # warp diffusion measure to template space first, then to MNI space
        antsApplyTransforms[
            '-d', '3',
            '-i', moving,
            '-o', output,
            '-r', mniTmp,
            '-t', warp2mni, trans2mni, warp2tmp, trans2tmp
        ] & FG
    
    return pjoin(templatePath, prefix + f'_InMNI_FA.nii.gz')



def sub2tmp2mni(templatePath, siteName, faImgs, N_proc):

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
    skel_mask= (skel.get_fdata()>0)*1.

    meanAttr=[]
    for faImg in faImgs:
        data= load(faImg).get_fdata()
        temp= data*skel_mask
        meanAttr.append(temp[temp>0].mean())

    return meanAttr


def main():

    parser = argparse.ArgumentParser(description='''Warps diffusion measures (FA, MD, GFA) to template space 
    and then to MNI space. Finally, calculates mean FA over IITmean_FA_skeleton.nii.gz''')
    parser.add_argument('-i', '--input', type=str, required=True,
        help='a .txt/.csv file that you used in/obtained from harmonization.py having two columns for (img,mask) pair. '
             'See documentation for more details')
    parser.add_argument('-s', '--site', type= str, required=True,
                        help='site name for locating template FA and mask in template directory')
    parser.add_argument('-t', '--template', type=str, required=True,
                        help='template directory where Mean_{site}_FA.nii.gz and {site}_Mask.nii.gz is located')
    parser.add_argument('--ncpu', help='number of cpus to use', default= '4')

    args = parser.parse_args()
    imgList=abspath(args.input)
    siteName=args.site
    templatePath=abspath(args.template)
    N_proc= int(args.ncpu)

    # read FA image list
    try:
        imgs, _ = read_imgs_masks(imgList)
        print('(Img,Mask) list is provided. FA images are assumed to be directoryOfImg/dti/ImgPrefix_FA.nii.gz, make sure they are there\n')
        faImgs= []

        for imgPath in imgs:
            directory = dirname(imgPath)
            prefix = basename(imgPath).split('.nii')[0]
            faImg= pjoin(directory, 'dti', prefix+ '_FA.nii.gz')
            if not isfile(faImg):
                raise FileNotFoundError(f'{faImg} not found. Did you run \"--create --debug\" and \"--process --debug\" before?')

            faImgs.append(faImg)


    except:
        faImgs= read_imgs(imgList)
        print('FA image list is provided.')


    # register and obtain *_InMNI_FA.nii.gz
    mniFAimgs= sub2tmp2mni(templatePath, siteName, faImgs, N_proc)
    
    # target harmonized
    if imgList.endswith('.modified.harmonized'):
        header = siteName+'_after'
    # reference
    elif imgList.endswith('.modified'):
        header = siteName
    # target unprocessed
    else:
        header = siteName+'_before'
        
    # FIXME user FA image list will use the header {siteName+'_before'}, which is not correct all the time
    # as shown in the above block:
    # reference should use {siteName} while harmonized target should use {siteName+'_after'}
    # impact of this discrepancy is minor since we deprecated use of FA image list
    
    outPrefix = pjoin(templatePath, header)
    
    print('\n\nComputing statistics\n\n')
    print(f'{siteName} site: ')
    site_means= analyzeStat(mniFAimgs)
    generate_csv(faImgs, site_means, outPrefix)
    
    # save statistics for future
    statFile= pjoin(templatePath, 'meanFAstat.csv')
    with open(statFile, 'a') as f:
        f.write(datetime.now().strftime('%m/%d/%y %H:%M')+',mean meanFA,std meanFA\n')
        f.write(f'{header},{np.mean(site_means)},{np.std(site_means)}\n')
        # print an empty line so future results, if appended, are visually separate
        f.write('\n')
            
    
    # print statistics on console
    print('\n\nPrinting statistics\n\n')
    with open(statFile) as f:
        print(f.read())

    # generate demonstrative plots
    ebar = harm_plot([site_means], [header], outPrefix)

    print(f'\nDetailed statistics, summary results, and demonstrative plots are saved in:\n\n{outPrefix}_stat.csv'
          f'\n{statFile}\n{ebar}\n')


if __name__ == '__main__':
    main()

