#!/usr/bin/env python

from os.path import abspath, dirname, join as pjoin
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
    # append 'lib/' so _version.py is discoverable
    FILEDIR = abspath(dirname(__file__))
    LIBDIR = dirname(FILEDIR)
    sys.path.append(LIBDIR)
    # get version info
    from _version import __version__

    # download test data
    test_data= sys.argv[1]
    test_unzip_dir= test_data.split('.')[0]

    chdir(pjoin(LIBDIR, 'tests'))
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



