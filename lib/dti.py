from plumbum.cmd import dtifit
from plumbum import FG

from util import *

def dti(imgPath, maskPath, inPrefix, outPrefix):

    print('dtifit ', imgPath)
    dtifit['-k', imgPath,
           '-m', maskPath,
           '-r', inPrefix+'.bvec',
           '-b', inPrefix+'.bval',
           '-o', outPrefix
            ] & FG

    vol= load(imgPath)
    mask= load(maskPath)
    masked_vol= applymask(vol.get_data(), mask.get_data())
    gfa_vol= gfa(masked_vol)
    save_nifti(outPrefix+'_GFA.nii.gz', gfa_vol, vol.affine, vol.header)


if __name__ == '__main__':
    pass