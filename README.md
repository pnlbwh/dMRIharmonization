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


# TODO
structuring element with radius 20 in cleanOutliers.py

# Nudge
Line 733 of `/home/tb571/anaconda3/lib/python3.6/site-packages/dipy/reconst/shm.py` edited
with `out[out>1.]= 1

# MATLAB call

rishFeatures_file('/home/tb571/Downloads/Harmonization-Python/connectom_prisma_demoData/A/connectom/dwi_A_connectom_st_b1200.nii.gz', '/home/tb571/Downloads/Harmonization-Python/connectom_prisma_demoData/A/connectom/mask.nii.gz', '/home/tb571/Downloads/Harmonization-Python/connectom_prisma_demoData/mat_lab/dwi_')



```
qb_model = QballModel(gtab, sh_order=N_shm, assume_normed= True, min_signal= eps)
    extend_bse= np.repeat(np.expand_dims(b0, 3), data.shape[3], axis=3)
    data=data/(extend_bse+eps)
    data[data<eps]=0
    data[data>1]=0

```

```
/home/tb571/Downloads/Harmonization-Python/lib/diff_map.py /home/tb571/Downloads/Harmonization-Python/connectom_prisma_demoData/py_thon/ /home/tb571/Downloads/Harmonization-Python/connectom_prisma_demoData/mat_lab/ 8

```
