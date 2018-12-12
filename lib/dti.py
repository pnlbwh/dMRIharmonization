from plumbum.cmd import dtifit
from plumbum import FG

import warnings
with warnings.catch_warnings():
    warnings.filterwarnings("ignore", category=FutureWarning)

    from dipy.io.image import load_nifti, save_nifti
    from dipy.reconst.odf import gfa
    from dipy.segment.mask import applymask

def dti(imgPath, maskPath, inPrefix, outPrefix):

    dtifit['-k', imgPath,
           '-m', maskPath,
           '-r', inPrefix+'.bvec',
           '-b', inPrefix+'.bval',
           '-o', outPrefix
            ] & FG

    vol, affine= load_nifti(imgPath)
    mask, _= load_nifti(maskPath)
    masked_vol= applymask(vol, mask)
    gfa_vol= gfa(masked_vol)
    save_nifti(outPrefix+'_GFA.nii.gz', gfa_vol, affine)


def main():
    dti('/home/tb571/Downloads/Harmonization-Python/connectom_prisma_demoData/A/connectom/dwi_A_connectom_st_b1200.nii.gz',
        '/home/tb571/Downloads/Harmonization-Python/connectom_prisma_demoData/A/connectom/mask.nii.gz',
        '/home/tb571/Downloads/Harmonization-Python/connectom_prisma_demoData/fa_test/abc')

if __name__ == '__main__':
    main()