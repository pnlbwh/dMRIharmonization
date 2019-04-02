![](doc/pnl-bwh-hms.png)

[![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.2584271.svg)](https://doi.org/10.5281/zenodo.2584271) [![Python](https://img.shields.io/badge/Python-3.6-green.svg)]() [![Platform](https://img.shields.io/badge/Platform-linux--64%20%7C%20osx--64-orange.svg)]()

*dMRIharmonization* repository is developed by Tashrif Billah, Sylvain Bouix, Suheyla Cetin Karayumak, and Yogesh Rathi, Brigham and Women's Hospital (Harvard Medical School).

Table of Contents
=================

   * [Table of Contents](#table-of-contents)
   * [Multi-site dMRI harmonization](#multi-site-dmri-harmonization)
   * [Citation](#citation)
   * [Dependencies](#dependencies)
   * [Installation](#installation)
      * [1. Install prerequisites](#1-install-prerequisites)
         * [Check system architecture](#check-system-architecture)
         * [Python 3](#python-3)
         * [FSL](#fsl)
         * [unringing](#unringing)
      * [2. Install pipeline](#2-install-pipeline)
         * [3. Configure your environment](#3-configure-your-environment)
   * [Running](#running)
   * [Tests](#tests)
      * [1. pipeline](#1-pipeline)
      * [2. unittest](#2-unittest)
   * [List of images](#list-of-images)
   * [Site names](#site-names)
   * [Multi threading](#multi-threading)
   * [Order of spherical harmonics](#order-of-spherical-harmonics)
   * [NRRD support](#nrrd-support)
   * [Preprocessing](#preprocessing)
         * [1. Denoising](#1-denoising)
         * [2. Bvalue mapping](#2-bvalue-mapping)
         * [3. Resampling](#3-resampling)
   * [Config](#config)
   * [Template](#template)
   * [List of outputs](#list-of-outputs)
      * [1. Folders](#1-folders)
      * [2. Files](#2-files)
   * [Template creation](#template-creation)
   * [Data harmonization](#data-harmonization)
   * [Debugging](#debugging)
         * [1. With the pipeline](#1-with-the-pipeline)
         * [2. Use seperately](#2-use-seperately)
   * [Travel heads](#travel-heads)
   * [Caveats](#caveats)
   * [Reference](#reference)

Created by [gh-md-toc](https://github.com/ekalinin/github-markdown-toc)


# Multi-site dMRI harmonization

Integrated study of multi-site diffusion MRI (dMRI) data can enable diagnosis, monitoring, and treatment of 
several brain diseases. However, data acquired on a variety of scanners cannot be integrated due to 
differences in acquisition parameters and scanner artifacts. Therefore, dMRI data has to be harmonized for 
joint analysis by removing scanner-specific differences. 
 
*dMRIharmonization* is a Python command line module that implements a method capable of removing scanner-specific effects. 
The proposed method eliminates inter-site variability in acquisition parameters, 
while preserves inter-subject anatomical variability. 

![](doc/flowchart.png)
See [Reference](#reference) for more details

The method leverages on rotation invariant spherical harmonics (RISH) features derived from spherical harmonic coefficients. 
In brief, the method constructs a scale map for each pair of reference and target sites. Applying the scale map to the 
target site, it aims to remove inter-site variability. 


# Citation

If this repository is useful in your research, please cite as below: 

Tashrif Billah, Sylvain Bouix, Suheyla Cetin Karayumak, and Yogesh Rathi *Multi-site Diffusion MRI Harmonization*, 
https://github.com/pnlbwh/dMRIharmoniziation, 2018, DOI: 10.5281/zenodo.2584271


# Dependencies

* ANTs
* FSL
* reisert/unring
* numpy
* scipy
* scikit-image
* dipy
* nibabel
* pynrrd


# Installation

## 1. Install prerequisites

Python 3, and FSL (ignore the one(s) you have already):

### Check system architecture

    uname -a # check if 32 or 64 bit

### Python 3

Download [Miniconda Python 3.6 bash installer](https://conda.io/miniconda.html) (32/64-bit based on your environment):
    
    sh Miniconda3-latest-Linux-x86_64.sh -b # -b flag is for license agreement

Activate the conda environment:

    source ~/miniconda3/bin/activate # should introduce '(base)' in front of each line

### FSL

Follow the [instruction](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FslInstallation) to download and install FSL.


### unringing

    git clone https://bitbucket.org/reisert/unring.git
    cd unring/fsl
    export PATH=$PATH:`pwd`
    unring --help

## 2. Install pipeline

Now that you have installed the prerequisite software, you are ready to install the pipeline:

    git clone https://github.com/pnlbwh/dMRIharmonization.git && cd dMRIharmonization
    conda create -f environmnet.yml     # you may comment out any existing package from environment.yml
    conda activate harmonization        # should introduce '(harmonization)' in front of each line


### 3. Configure your environment

Make sure the following executables are in your path:

    antsMultivariateTemplateConstruction2.sh
    antsApplyTransforms
    antsRegistrationSyNQuick.sh
    dtifit
    unring.a64
    
You can check them as follows:

    which dtifit
    
If any of them does not exist, add that to your path:

    export PATH=$PATH:/directory/of/executable
    
`conda` should already put the ANTs scripts in your path. However, if you choose to use pre-installed ANTs scripts, 
you may need to define [ANTSPATH](https://github.com/ANTsX/ANTs/wiki/Compiling-ANTs-on-Linux-and-Mac-OS#set-path-and-antspath)


# Running

Upon successful installation, you should be able to see the help message

`$ lib/harmonization.py --help`

    Template creation, harmonization, and debugging
    
    Usage:
        harmonization.py [SWITCHES] 
    
    Meta-switches:
        -h, --help                          Prints this help message and quits
        --help-all                          Prints help messages of all sub-commands and quits
        -v, --version                       Prints the program's version and quits
    
    Switches:
        --bvalMap VALUE:str                 specify a bmax to scale bvalues into
        --create                            turn on this flag to create template
        --debug                             turn on this flag to debug harmonized data (valid only with --process)
        --denoise                           turn on this flag to denoise voxel data
        --force                             turn on this flag to overwrite existing data
        --harm_list VALUE:ExistingFile      harmonized csv/txt file with first column for dwi and 2nd column for mask: dwi1,mask1
                                            dwi2,mask2 ...
        --nproc VALUE:str                   number of processes/threads to use (-1 for all available, may slow down your system);
                                            the default is 8
        --nshm VALUE:str                    spherical harmonic order; the default is 6
        --process                           turn on this flag to harmonize
        --ref_list VALUE:ExistingFile       reference csv/txt file with first column for dwi and 2nd column for mask: dwi1,mask1
                                            dwi2,mask2 ...
        --ref_name VALUE:str                reference site name; required
        --resample VALUE:str                voxel size MxNxO to resample into
        --tar_list VALUE:ExistingFile       target csv/txt file with first column for dwi and 2nd column for mask: dwi1,mask1
                                            dwi2,mask2 ...; required
        --tar_name VALUE:str                target site name; required
        --template VALUE:str                template directory; required
        --travelHeads                       travelling heads


The `harmonization.py` cli takes in the following arguments that are explained below.

* [--bvalMap VALUE:str](#2-bvalue-mapping)
* [--create](#template-creation)
* [--debug](#debugging)
* [--denoise](#1-denoising)
* --force
* [--harm_list VALUE:ExistingFile](#2.-use-seperately)
* [--nproc VALUE:str](#multi-threading)
* [--nshm VALUE:str](#order-of-spherical-harmonics)
* [--process](#data-harmonization)
* [--ref_list VALUE:ExistingFile](#list-of-images)
* [--ref_name VALUE:str](#site-names)
* [--resample VALUE:str](#3-resampling)
* [--tar_list VALUE:ExistingFile](#list-of-images)
* [--tar_name VALUE:str](#site-names)
* [--template VALUE:str](#template)
* [--travelHeads](#travel-heads)



# Tests

A small test data is provided with each [release](https://github.com/pnlbwh/Harmonization-Python/releases). 

## 1. pipeline
You may test the whole pipeline as follows:
    
    cd lib/tests
    ./test_pipeline.sh
    
NOTE: running the above test should take roughly an hour.

`./pipeline_test.sh` will download test data, and run the whole processing pipeline on them. 
If the test is successful and complete, you should see the following output on the command line. 
    
    CONNECTOM mean FA:  0.9115503529641469
    PRISMA mean FA before harmonization:  0.8011634268930871
    PRISMA mean FA after harmonization:  0.8222395983556193

## 2. unittest
You may run smaller and faster unittest as follows:
    
    python -m unittest discover -v lib/tests/    
    

# List of images

Two image lists (.txt or .csv) are required for the reference and target sites:

    --ref_list cidar_imgs_masks.csv
    --tar_list mrc_imgs_masks.txt
    
Each `\n` delimited line in the list file should have `,` seperated dwi and mask:
 
    case001_dwi.nii.gz,case001_mask.nii.gz
    case002_dwi.nrrd,case002_mask.nii.gz
    case003_dwi.nii,case003_mask.nhdr
    ...
    ...

It is assumed that *.bval* and *.bvec* are in the same location 
and have the same prefix as that of the `dwi.nii.gz` image.

During template construction, you want to provide a **small list of images** to `--ref_list` and `--tar_list`. However, 
provide a **complete list of images** you want to harmonize to `--tar_list` during [data harmonization](#data-harmonization)

# Site names

Site names are used to properly identify files in the template:
    
     --ref_name CIDAR
     --tar_name MRC
     

# Multi threading

Processing can be multi-threaded over the cases. Besides, `antsMultivariateTemplateConstruction2.sh` utilizes 
multiple threads to speed-up template construction. 

    --nproc 8 # default is 4, use -1 for all available
   
However, multi-threading comes with a price of slowing down other processes that may be running in your system. So, it 
is advisable to leave out at least two cores for other processes to run smoothly.


# Order of spherical harmonics

RISH features are derived from spherical harmonic coefficients. The order of spherical harmonic coefficients you can use 
is limited by the lowest number of gradients present in the diffusion images between reference and target sites. 


| `--nshm` | # of shm coefficients | Required # gradients >= |
|----------|-----------------------|-------------------------|
|    0     |            1          |           1             |
|    2     |            6          |           6             |
|    4     |            15         |           15            |
|    6     |            28         |           28            |
|    8     |            45         |           45            |


# NRRD support

The pipeline is written for NIFTI image format. However, NRRD support is incorporated through [NIFTI --> NRRD](https://github.com/pnlbwh/Harmonization-Python/blob/parallel/lib/preprocess.py#L78) 
conversion on the fly.

See Tashrif Billah, Isaiah Norton, Ryan Eckbo, and Sylvain Bouix, Processing pipeline for anatomical and diffusion weighted images, 
https://github.com/pnlbwh/pnlpipe, 2018, DOI: 10.5281/zenodo.2584271 for more details on the conversion method.


# Preprocessing

There are three preprocessing steps. More than one step can be specified at a time. However, they are performed in 
the following order.

### 1. Denoising
    
    --denoise        # turn on this flag to denoise voxel data

### 2. Bvalue mapping

    --bvalMap VALUE  # specify a bmax to scale bvalues into    

### 3. Resampling

    --resample VALUE # voxel size MxNxO to resample into


After preprocessing, the image lists are saved with `.modified` extension in the same location of provided lists, 
and used for further processing.
 
# Config   
    
With all the arguments passed to `harmonization.py` cli, a *config.ini* file is created in `lib/`. This file is used by 
submodules of the harmonization pipeline during executation. A sample `config.ini` file looks like the following:

    [DEFAULT]
    N_shm = 6
    N_proc = 8
    resample = 1.5x1.5x1.5
    bvalMap = 1000
    denoise = 0
    travelHeads = 0
    debug = 1
    diffusionMeasures = MD,FA,GFA

# Template

Template directory is passed by `--template` arguments. It is primarily used by `antsMultivariateTemplateConstruction2.sh` 
to save transform files. Then, various scales maps (*Scale_\*.nii.gz*) and mean templates (*Mean_\*.nii.gz*) 
are created in this directory. There are other template files: 
*Delta_\*.nii.gz*, *PercentageDiff_\*.nii.gz*, and *PercentageDiff_\*smooth.nii.gz*. You may look at them to know how good 
has been the created template.
    
# List of outputs

Several files are created down the pipeline. They are organized with proper folder hierarchy and naming:
    
## 1. Folders

In each input image directory, two folders are created: `dti` and `rish`. The `dti` folder stores diffusion measures 
and corresponding transform files. On the other hand `rish` folder stores RISH features and corresponding transform files.

## 2. Files
    
The harmonization pipeline will generate a bunch of diffusion measures and RISH features. All 
the measures are saved as *.nii.gz* format with associated *.bval* and *.bvec* where necessary. 

Firstly, see the files with prefix *harmonized_* in target site image directories. They are the harmonized 
diffusion weighted images.

Secondly, you can look at files with suffix *_denoised.nii.gz*, *_bmapped.nii.gz*, and *_resampled.nii.gz*. 
Denoising, bvalue mapping, and resampling are performed in the above order and files are saved with the last suffix 
unless `--debug` option is used. In the latter case, files after each preprocessing step are saved.

Finally, there are other files:
    
|         string         |             description                      |
|------------------------|----------------------------------------------|
|    \*Warped\*.nii.gz   |  warped to template space                    |
|    \*_InMNI\*.nii.gz   |  warped to MNI space                         |
|    ToSubjectSpace_\*   |  registration of template to subject space   |

      
# Template creation
    
RISH features are diffusion measures are computed from all the input images. Images of two modalites: 
*FA* and *L0* (RISH order 0) are provided as input to `antsMultivariateTemplateConstruction2.sh`. A `template0.nii.gz` is 
created. Afterwards, various scales maps (*Scale_\*.nii.gz*) and mean templates (*Mean_\*.nii.gz*) 
are created in `--template` directory. See [Reference](#referece) for details on the template construction method 
and [Template](#Template) for list of outputs. A sample template construction command is given below:

    lib/harmonization.py \
    --ref_list test_data/ref_caselist.txt \
    --tar_list test_data/target_caselist.txt \
    --ref_name REF \
    --tar_name TAR \
    --template test_data/template/ \
    --nshm 6 \
    --bvalMap 1000 \
    --resample 1.5x1.5x1.5 \
    --nproc 8 \
    --create

Note: Replace each of the above paths with absolute path as needed.

During template construction, you want to provide a **small list of images** to `--ref_list` and `--tar_list`. However, 
provide a **complete list of images** you want to harmonize to `--tar_list` during [data harmonization](#data-harmonization)


# Data harmonization

After template creation is completed, there should be scales maps (*Scale_\*.nii.gz*) in the template directory. 
The scale maps are used to scale RISH features of the reference site and then reconstruct diffusion weighted images. 
See [Reference](#referece) for details on the signal reconstruction.
    
    lib/harmonization.py \
    --tar_list test_data/target_caselist.txt \
    --tar_name TAR \
    --template test_data/template/ \
    --nshm 6 \
    --bvalMap 1000 \
    --resample 1.5x1.5x1.5 \
    --nproc 8 \
    --process

Note: Replace each of the above paths with absolute path as needed.

During template construction, you want to provide a **small list of images** to `--ref_list` and `--tar_list`. However, 
provide a **complete list of images** you want to harmonize to `--tar_list` during [data harmonization](#data-harmonization)


# Debugging

dMRIharmonization toolbox is provided with a debugging capability to test how good has been the 
harmonization. `--debug` can be run with any (or all) of `--create` and `--process` options or separately. 

Once harmonization is done, there are three types of data:

* reference site 
* target site
* harmonized target

If the three data sets are brought into a common space, then looking at the mean FA over the whole brain, 
we can comment on the goodness of harmonization. If data is properly harmonized, then mean FA of the harmonized target 
should be almost equal to the mean FA of that of reference site.

In details:

* reference data is proprocessed, and registered to reference template space and then to MNI space ([IITmean_FA.nii.gz](http://www.iit.edu/~mri/IITHumanBrainAtlas.html))
* unprocessed target data is directly registered to MNI space
* harmonized target data is registered to target template space and then to MNI space
* once the data are in MNI space, we calculate mean FA over the [IITmean_FA_skeleton.nii.gz](http://www.iit.edu/~mri/IITHumanBrainAtlas.html) 
 

The numbers should like like below:

    Printing statistics :
    CIDAR mean FA:  0.5217237675408243
    BSNIP mean FA before harmonization:  0.5072286796848892
    BSNIP mean FA after harmonization:  0.5321998242139347

As we see above, BSNIP (target) mean FA over [IITmean_FA_skeleton.nii.gz](http://www.iit.edu/~mri/IITHumanBrainAtlas.html) 
for all the site images after harmonization increased to be almost equal to that of CIDAR (reference) mean FA.

Now there are two ways to debug:

### 1. With the pipeline 
Use `--debug` flag with any (or all) of `--create` and `--process`

### 2. Use seperately 
If you would like to debug at a later time, you need to specify three images lists:

* `--ref_list`: use the reference list with `.modified` extension
* `--tar_list`: use the unprocessed target list **without** the `.modified` extension
* `--harm_list`: use the harmonized target list with `.harmonized` extension

    lib/harmonization.py --ref_list ref.txt.modified --tar_list target.csv --harm_list target.csv.modified.harmonized

NOTE: You should run the pipeline first before debugging separately because `--debug` makes use of files created 
in the pipeline.

# Travel heads

Traveling heads mean *same set of subjects scanned with different scanners*. 
For example, scanning same set of people in Brigham and Women's Hospital and 
Massachusetts General Hospital. Therefore, users should set `--travelHeads` 
if they have same set of subjects across sites.

# Caveats

`antsMultivariateTemplateConstruction2.sh`: all the images need to have unique
prefix because transform files are created in the same `--template ./template/` directory. The value of `--template` 
should have `/` at the end. The pipeline appends one if there is not, but it is good to include it when specifying.

# Reference

Karayumak, S.C., Bouix, S., Ning, L., James, A., Crow, T., Shenton, M., Kubicki, M. and Rathi, Y., 2019. 
Retrospective harmonization of multi-site diffusion MRI data acquired with different acquisition parameters. 
Neuroimage, 184, pp.180-200.