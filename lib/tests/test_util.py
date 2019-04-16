from os.path import abspath, dirname, basename, exists, join as pjoin, isfile, split as psplit
from os import mkdir
import sys
from shutil import rmtree, copyfile
import unittest
from subprocess import check_call

FILEDIR= abspath(dirname(__file__))
LIBDIR= dirname(FILEDIR)

# sys.path.append(FILEDIR)
sys.path.append(LIBDIR)

import warnings
with warnings.catch_warnings():
    warnings.filterwarnings("ignore", category=FutureWarning)
    from nibabel import load, Nifti1Image
    from dipy.io import read_bvals_bvecs

def save_nifti(fname, data, affine, hdr=None):
    if data.dtype.name=='uint8':
        hdr.set_data_dtype('uint8')
    else:
        hdr.set_data_dtype('float32')

    result_img = Nifti1Image(data, affine, header=hdr)
    result_img.to_filename(fname)

from numpy import array
from conversion import write_bvals

test_data= 'connectom_prisma.zip'
test_unzip_dir= pjoin(FILEDIR, test_data.split('.')[0])
if not exists(pjoin(FILEDIR, test_unzip_dir)):
    download_script= pjoin(FILEDIR, 'download_data.py')
    with open(download_script) as f:
        check_call(f'{download_script} {test_data}', shell= True)