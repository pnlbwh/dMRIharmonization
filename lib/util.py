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

import os, shutil, configparser, sys
import numpy as np
from subprocess import check_call

import warnings
with warnings.catch_warnings():
    warnings.filterwarnings("ignore", category=FutureWarning)
    from nibabel import load, Nifti1Image
    from dipy.io.image import load_nifti
    from dipy.io import read_bvals_bvecs
    from dipy.core.gradients import gradient_table
    from dipy.reconst.shm import QballModel
    from dipy.reconst.odf import gfa
    from dipy.segment.mask import applymask
    import dipy.reconst.dti as dipyDti


def save_nifti(fname, data, affine, hdr=None):
    if data.dtype.name=='uint8':
        hdr.set_data_dtype('uint8')
    else:
        hdr.set_data_dtype('float32')

    result_img = Nifti1Image(data, affine, header=hdr)
    result_img.to_filename(fname)


def convertedPath(imgPath):

    if imgPath.endswith('.nhdr') or imgPath.endswith('.nrrd'):
        imgPath = imgPath.split('.')[0] + '.nii.gz'

    return imgPath


B0_THRESH= 50.
B_QUANT= 50.
BSHELL_MIN_DIST= 100.
