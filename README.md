# Harmonization-Python
Refactoring MATLAB Harmonization tool into Python


# Template creation

'''
python -m pdb \
/home/tb571/Downloads/Harmonization-Python/lib/harmonization.py \
--reference /home/tb571/Downloads/Harmonization-Python/test_data/ref_caselist.txt \
--target /home/tb571/Downloads/Harmonization-Python/test_data/target_caselist.txt \
--template /home/tb571/Downloads/Harmonization-Python/test_data/template/ \
--travelHeads \
--force \
--N_shm 6 \
--create

'''


# Harmonization

'''
python -m pdb \
/home/tb571/Downloads/Harmonization-Python/lib/harmonization.py \
--N_shm 6 \
--denoise \
--force \
--process \
--target /home/tb571/Downloads/Harmonization-Python/test_data/target_caselist.txt \
--template /home/tb571/Downloads/Harmonization-Python/test_data/template/ \
--travelHeads
'''


artifact of antsMultivariateTemplateConstruction2.sh: reference images need to have unique
prefix. See lines 42 and 43 in buildTemplate.py

--template /home/tb571/Downloads/Harmonization-Python/test_data/template/: the last / is 
important, it is an artifact of antsMultivariateTemplateConstruction2.sh
