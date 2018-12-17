#!/usr/bin/bash

# set -x

sites=(connectom prisma)
cases=(a b c)

for s in ${sites[@]}
do
    for c in ${cases[@]}
    do
        echo "$s, $c"
        /home/tb571/Downloads/Harmonization-Python/lib/test_data_creation.py \
        --prefix /home/tb571/Downloads/Harmonization-Python/test_data/test_$c/$s/${s}_${c}_dwi \
        --nifti /home/tb571/Downloads/Harmonization-Python/connectom_prisma_demoData/${c^^}/$s/"dwi_${c^^}_${s}_st_b1200.nii.gz" \
        --bval /home/tb571/Downloads/Harmonization-Python/connectom_prisma_demoData/${c^^}/$s/"dwi_${c^^}_${s}_st_b1200.bval" \
        --bvec /home/tb571/Downloads/Harmonization-Python/connectom_prisma_demoData/${c^^}/$s/"dwi_${c^^}_${s}_st_b1200.bvec"
        bet /home/tb571/Downloads/Harmonization-Python/test_data/test_$c/$s/"${s}_${c}_dwi.nii.gz" \
        /home/tb571/Downloads/Harmonization-Python/test_data/test_$c/$s/"${s}_${c}_mask.nii.gz"
    done
done

# set +x