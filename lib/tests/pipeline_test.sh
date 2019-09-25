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
    wget https://github.com/pnlbwh/Harmonization-Python/releases/download/v${__version__}/${test_data}.zip
fi

tar -xzvf ${test_data}.zip

cd ${test_data}
CURRDIR=`pwd`


# append path to image list and write back
write_list connectom.txt
write_list prisma.txt


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
--create --process --debug

