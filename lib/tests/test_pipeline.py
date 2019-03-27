#!/usr/bin/env python

import os, sys, shutil
from subprocess import check_call

# append path so _version.py is discoverable
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
# get version info
from _version import __version__

# download test data
test_data= 'connectom_prisma.zip'
# check_call(['wget', f'https://github.com/pnlbwh/Harmonization-Python/releases/download/v{__version__}/{test_data}'])
# check_call(' '.join(['tar', '-xzvf', f'{test_data}']), shell= True)

os.chdir('connectom_prisma')
CURRDIR= os.getcwd()


def write_list(caselist):

    outlist= '/tmp/harm_caselist.txt'
    with open(outlist, 'w') as fw:
        with open(caselist) as f:
            content = f.read()
            for line, row in enumerate(content.split()):
                temp = [element for element in row.split(',') if element]  # handling w/space
                fw.write(f'{CURRDIR}/{temp[0]},{CURRDIR}/{temp[1]}\n')

    shutil.move(outlist, caselist)



# # append path to image list and write back
write_list('connectom.txt')
write_list('prisma.txt')


# run test
check_call((' ').join(['../../harmonization.py',
                      '--bvalMap', '1000',
                      '--resample', '1.5x1.5x1.5',
                      '--template', './template/',
                      '--ref_list', 'connectom.txt',
                      '--tar_list', 'prisma.txt',
                      '--ref_name', 'CONNECTOM',
                      '--tar_name', 'PRISMA',
                      '--nshm', '4',
                      '--nproc', '-1',
                      '--create', '--process', '--debug']), shell= True)