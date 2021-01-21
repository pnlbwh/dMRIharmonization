#!/usr/bin/env bash

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


EXIT()
{
    echo ''
    echo $1
    exit 1
}    


write_list()
{
    CASELIST=$1
    TMPFILE=/tmp/harm_caselist.txt
    while IFS=, read -r img mask; do echo $CURRDIR/$img,$CURRDIR/$mask; done < $1 > $TMPFILE

    mv $TMPFILE $CASELIST
}

# get version info
IFS=" = ", read -r _ v < ../../_version.py
__version__=`echo $v | xargs`

# download test data
test_data=connectom_prisma # change this value if test data name is changed
if [ ! -f ${test_data}.zip ]
then
    wget https://github.com/pnlbwh/dMRIharmonization/releases/download/v${__version__}/${test_data}.zip
fi

tar -xzvf ${test_data}.zip

cd ${test_data}
CURRDIR=`pwd`


# append path to image list and write back
write_list connectom.txt
write_list prisma.txt


### Run pipeline and obtain statistics when same number of matched reference and target images are used in
### tempalate creation and harmonization
# run test
../../harmonization.py \
--bvalMap 1000 \
--resample 1.5x1.5x1.5 \
--template ./template/ \
--ref_list connectom.txt \
--tar_list prisma.txt \
--ref_name CONNECTOM \
--tar_name PRISMA \
--nproc -1 \
--create --process --debug || EXIT 'harmonization.py with --create --process --debug failed'
exit
# recompute statistics
../../harmonization.py \
--template ./template/ \
--ref_list connectom.txt.modified \
--tar_list prisma.txt \
--harm_list prisma.txt.modified.harmonized \
--ref_name CONNECTOM \
--tar_name PRISMA \
--stats || EXIT 'harmonization.py with --stats failed'
# ===============================================================================================================



### Run pipeline and obtain statistics when small set of matched reference and target images are used in template creation
### and a larger set (does not have to be mutually exclusive from the former) of target images are used in harmonization

# test the following advanced parameter
export TEMPLATE_CONSTRUCT_CORES=6

# --create and --debug block
../../harmonization.py \
--bvalMap 1000 \
--resample 1.5x1.5x1.5 \
--template ./template/ \
--ref_list connectom.txt \
--tar_list prisma.txt \
--ref_name CONNECTOM \
--tar_name PRISMA \
--travelHeads \
--nproc -1 \
--create --debug --force || EXIT 'harmonization.py with --create --debug --force failed'

../../harmonization.py \
--bvalMap 1000 \
--resample 1.5x1.5x1.5 \
--template ./template/ \
--ref_list connectom.txt \
--tar_list prisma.txt \
--ref_name CONNECTOM \
--tar_name PRISMA \
--travelHeads \
--nproc -1 \
--create --debug || EXIT 'harmonization.py with --create --debug failed'

# --process and --debug block
../../harmonization.py \
--bvalMap 1000 \
--resample 1.5x1.5x1.5 \
--template ./template/ \
--tar_list prisma.txt \
--tar_name PRISMA \
--ref_list connectom.txt \
--nproc -1 \
--process --debug --force || EXIT 'harmonization.py with --process --debug --force failed'

../../harmonization.py \
--bvalMap 1000 \
--resample 1.5x1.5x1.5 \
--template ./template/ \
--tar_list prisma.txt \
--tar_name PRISMA \
--ref_list connectom.txt \
--nproc -1 \
--process --debug || EXIT 'harmonization.py with --process --debug failed'

# ===============================================================================================================

# same bvalue, resolution block
cp connectom.txt.modified connectom_same.txt
cp prisma.txt.modified prisma_same.txt
../../harmonization.py \
--template ./template/ \
--ref_list connectom_same.txt \
--tar_list prisma_same.txt \
--ref_name CONNECTOM \
--tar_name PRISMA \
--nproc -1 \
--create --process || EXIT 'harmonization.py for same bvalue, resolution with --create --process failed'


# ===============================================================================================================
# compute statistics
../fa_skeleton_test.py -i connectom.txt.modified \
-s CONNECTOM -t template/ || EXIT 'fa_skeleton_test.py failed for modified reference'

../fa_skeleton_test.py -i prisma.txt \
-s PRISMA -t template/ || EXIT 'fa_skeleton_test.py failed for given target'

../fa_skeleton_test.py -i prisma.txt.modified.harmonized \
-s PRISMA -t template/ || EXIT 'fa_skeleton_test.py failed for harmonized target'



# now that all files are created, recompute statistics
../../harmonization.py \
--template ./template/ \
--ref_list connectom.txt.modified \
--tar_list prisma.txt \
--harm_list prisma.txt.modified.harmonized \
--ref_name CONNECTOM \
--tar_name PRISMA \
--stats || EXIT 'harmonization.py with --stats failed'

# ===============================================================================================================


