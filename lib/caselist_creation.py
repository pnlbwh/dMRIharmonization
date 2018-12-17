#!/usr/bin/env python
import argparse, glob, os

def main():

    parser = argparse.ArgumentParser(description='create caselist for test data')
    parser.add_argument('-d', '--dir', type=str, required=True,
                        help='directory of test_data')

    parser.add_argument('-s', '--subDir', type=str, required=False,
                    help='sub directory inside each case')

    parser.add_argument('-o', '--out', type=str, required=True,
                help='directory of test_data')

    args = parser.parse_args()
    directory= os.path.abspath(args.dir)
    subDir= args.subDir
    if not subDir:
        subDir= '.'
    output= args.out

    with open(output,'w') as f:
        cases= [d for d in os.listdir(directory) if os.path.isdir(d)]
        for c in cases:

            searchPath= os.path.join(directory, c, subDir, '*dwi*.nii.gz')
            imgFile= glob.glob(searchPath)
            if not len(imgFile):
                searchPath= os.path.join(directory, c, subDir, '*DWI*.nii.gz')
                imgFile= glob.glob(searchPath)
            if len(imgFile)>1:
                raise AttributeError(f'Multiple dwis exist in {searchPath}')


            searchPath= os.path.join(directory, c, subDir, '*mask*.nii.gz')
            maskFile= glob.glob(searchPath)
            if not len(maskFile):
                searchPath= os.path.join(directory, c, subDir, '*MASK*.nii.gz')
                maskFile= glob.glob(searchPath)
            if len(maskFile)>1:
                raise AttributeError(f'Multiple masks exist in {searchPath}')

            f.write(f'{imgFile[0]},{maskFile[0]}\n')


if __name__=='__main__':
    main()

'''
/home/tb571/Downloads/Harmonization-Python/lib/caselist_creation.py -d ./ -s connectom -o ./ref_caselist.txt
/home/tb571/Downloads/Harmonization-Python/lib/caselist_creation.py -d ./ -s prisma -o ./target_caselist.txt
'''