#!/usr/bin/env bash


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
wget https://github.com/pnlbwh/Harmonization-Python/releases/download/v${__version__}/${test_data}.zip
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
--nshm 4 \
--nproc -1 \
--create --process --debug