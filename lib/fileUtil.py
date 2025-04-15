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

from util import exists, warnings, makedirs, rmtree, pjoin, dirname

def check_csv(file, force):

    with open(file) as f:
        content= f.read()

        for line, row in enumerate(content.split()):
            dwi_mask= [element for element in row.split(',') if element] # handling w/space
            if len(dwi_mask) != 2:
                raise FileNotFoundError(f'Columns don\'t have same number of entries: check line {line} in {file}')

            dirCheckFlag= 1
            for img in dwi_mask:
                if not exists(img):
                    raise FileNotFoundError(f'{img} does not exist: check line {line} in {file}')

                elif dirCheckFlag:
                    # create DTI and harmonization directory
                    dtiPath= pjoin(dirname(img),'dti')
                    check_dir(dtiPath, force)

                    harmPath= pjoin(dirname(img),'harm')
                    check_dir(harmPath, force)

                    dirCheckFlag= 0


def check_dir(path, force):
    if exists(path) and force:
        warnings.warn(f'{path} exists and will be overwritten')
        rmtree(path)
        makedirs(path)
    elif not exists(path):
        makedirs(path)
    else:
        warnings.warn(f'{path} exists, --force not specified, continuing with existing directory')


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

# multi-shell-dMRIharmonization takes NIFTI input only
# this block may be uncommented in a future design
# convert NRRD to NIFTI on the fly
# from conversion import nifti_write
# def nrrd2nifti(imgPath):
#
#     if imgPath.endswith('.nrrd'):
#         niftiImgPrefix= imgPath.split('.nrrd')[0]
#     elif imgPath.endswith('.nhdr'):
#         niftiImgPrefix= imgPath.split('.nhdr')[0]
#     else:
#         return imgPath
#
#     nifti_write(imgPath, niftiImgPrefix)
#     return niftiImgPrefix+'.nii.gz'

