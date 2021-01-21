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

from os.path import abspath, dirname, join as pjoin, isfile
from os import chdir, getcwd
import sys, shutil
from subprocess import check_call

def write_list(caselist):

    outlist= '/tmp/harm_caselist.txt'
    with open(outlist, 'w') as fw:
        with open(caselist) as f:
            content = f.read()
            for line, row in enumerate(content.split()):
                temp = [element for element in row.split(',') if element]  # handling w/space
                fw.write(f'{CURRDIR}/{temp[0]},{CURRDIR}/{temp[1]}\n')

    shutil.move(outlist, caselist)


def main():

    FILEDIR = abspath(dirname(__file__))
    LIBDIR = dirname(FILEDIR)
    ROOTDIR = dirname(LIBDIR)
    # append 'lib/' so pipeline functions are discoverable
    sys.path.append(LIBDIR)
    # append 'dMRIharmonization/' so _version.py is discoverable
    sys.path.append(ROOTDIR)
    # get version info
    from _version import __version__

    if len(sys.argv)==1 or sys.argv[1] in ['-h', '--help']:
        print('Usage: download_data.py dataName.zip\n'
              f'dataName.zip is the name of the zip file in https://github.com/pnlbwh/Harmonization-Python/releases/latest')
        return

    # download test data
    test_data= sys.argv[1]
    test_unzip_dir= test_data.split('.')[0]

    chdir(pjoin(LIBDIR, 'tests'))
    if not isfile(test_data):
        check_call(['wget', f'https://github.com/pnlbwh/Harmonization-Python/releases/download/v{__version__}/{test_data}'])
    check_call(' '.join(['tar', '-xzvf', f'{test_data}']), shell=True)

    chdir(test_unzip_dir)
    global CURRDIR
    CURRDIR= getcwd()

    # append path to image list and write back
    write_list('connectom.txt')
    write_list('prisma.txt')

if __name__ == '__main__':
    main()



