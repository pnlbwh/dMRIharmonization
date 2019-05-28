![](doc/pnl-bwh-hms.png)

[![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.2584275.svg)](https://doi.org/10.5281/zenodo.2584275) [![Python](https://img.shields.io/badge/Python-3.6-green.svg)]() [![Platform](https://img.shields.io/badge/Platform-linux--64%20%7C%20osx--64-orange.svg)]()

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
         * [unringing](#unringing)
      * [2. Install pipeline](#2-install-pipeline)
      * [3. Download IIT templates](#3-download-iit-templates)
      * [4. Configure your environment](#4-configure-your-environment)
   * [Running](#running)
   * [Tests](#tests)
      * [1. pipeline](#1-pipeline)
      * [2. unittest](#2-unittest)
   * [List of images](#list-of-images)
   * [Site names](#site-names)
   * [Multi threading](#multi-threading)
   * [Order of spherical harmonics](#order-of-spherical-harmonics)
   * [Number of zero padding](#number-of-zero-padding)
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
      * [3. With a list of FA images](#3-with-a-list-of-fa-images)
   * [Travel heads](#travel-heads)
   * [Caveats/Issues](#caveatsissues)
      * [1. Template path](#1-template-path)
      * [2. Multi-processing](#2-multi-processing)
      * [3. Tracker](#3-tracker)
   * [Reference](#reference)



Table of Contents created by [gh-md-toc](https://github.com/ekalinin/github-markdown-toc)



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

If this repository is useful in your research, please cite all of the following: 

* Billah T*, Cetin Karayumak S*, Bouix S, Rathi Y. Multi-site Diffusion MRI Harmonization, 
https://github.com/pnlbwh/dMRIharmoniziation, 2019, doi: 10.5281/zenodo.2584275

    \* *denotes equal first authroship*


* Cetin Karayumak S, Bouix S, Ning L, James A, Crow T, Shenton M, Kubicki M, Rathi Y. Retrospective harmonization of multi-site diffusion MRI data 
acquired with different acquisition parameters. Neuroimage. 2019 Jan 1;184:180-200. 
doi: 10.1016/j.neuroimage.2018.08.073. Epub 2018 Sep 8. PubMed PMID: 30205206; PubMed Central PMCID: PMC6230479.


* Mirzaalian H, Ning L, Savadjiev P, Pasternak O, Bouix S, Michailovich O, Karmacharya S, Grant G, Marx CE, Morey RA, Flashman LA, George MS, 
McAllister TW, Andaluz N, Shutter L, Coimbra R, Zafonte RD, Coleman MJ, Kubicki M, Westin CF, Stein MB, Shenton ME, Rathi Y. 
Multi-site harmonization of diffusion MRI data in a registration framework. Brain Imaging Behav. 2018 Feb;12(1):284-295. 
doi:10.1007/s11682-016-9670-y. PubMed PMID: 28176263.



# Dependencies

* ANTs = 2.2.0
* reisert/unring
* FSL = 5.0.11
* numpy = 1.16.2
* scipy = 1.2.1
* scikit-image = 0.15.0
* dipy = 0.16.0
* nibabel = 2.3.0
* pynrrd = 0.3.6
* conversion = 2.0

**NOTE** The above versions were used for developing the repository. However, *dMRIharmonization* should work on 
any advanced version. 


# Installation

## 1. Install prerequisites

You may ignore installation instruction for any software module that you have already.

### Check system architecture

    uname -a # check if 32 or 64 bit

### Python 3

Download [Miniconda Python 3.6 bash installer](https://docs.conda.io/en/latest/miniconda.html) (32/64-bit based on your environment):
    
    sh Miniconda3-latest-Linux-x86_64.sh -b # -b flag is for license agreement

Activate the conda environment:

    source ~/miniconda3/bin/activate # should introduce '(base)' in front of each line

### unringing

    git clone https://bitbucket.org/reisert/unring.git
    cd unring/fsl
    export PATH=$PATH:`pwd`
    unring --help

## 2. Install pipeline

Now that you have installed the prerequisite software, you are ready to install the pipeline:

    git clone https://github.com/pnlbwh/dMRIharmonization.git && cd dMRIharmonization
    conda env create -f environment.yml    # you may comment out any existing package from environment.yml
    conda activate harmonization           # should introduce '(harmonization)' in front of each line


Alternatively, if you already have ANTs, you can continue using your python environment by directly installing 
the prerequisite libraries:

    pip install -r requirements.txt --upgrade


## 3. Download IIT templates

dMRIharmonization toolbox is provided with a debugging capability to test how good has been the 
harmonization. For debug to work and **tests** to run, download the following data from [IIT HUMAN BRAIN ATLAS](http://www.iit.edu/~mri/IITHumanBrainAtlas.html) 
and place them in `IITAtlas/` directory:

* [IITmean_FA.nii.gz](https://www.nitrc.org/frs/download.php/6898/IITmean_FA.nii.gz) 
* [IITmean_FA_skeleton.nii.gz](https://www.nitrc.org/frs/download.php/6897/IITmean_FA_skeleton.nii.gz)


## 4. Configure your environment

Make sure the following executables are in your path:

    antsMultivariateTemplateConstruction2.sh
    antsApplyTransforms
    antsRegistrationSyNQuick.sh
    unring.a64
    
You can check them as follows:

    which dtifit
    
If any of them does not exist, add that to your path:

    export PATH=$PATH:/directory/of/executable
    
`conda activate harmonization` should already put the ANTs scripts in your path. Yet, you should set ANTSPATH as follows:
    
    export ANTSPATH=~/miniconda3/envs/harmonization/bin/

However, if you choose to use pre-installed ANTs scripts, you can define ANTSPATH according to [this](https://github.com/ANTsX/ANTs/wiki/Compiling-ANTs-on-Linux-and-Mac-OS#set-path-and-antspath) instruction.


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
        --nzero VALUE:str                   number of zero padding for denoising skull region during signal reconstruction; the default is 10
        --process                           turn on this flag to harmonize
        --ref_list VALUE:ExistingFile       reference csv/txt file with first column for dwi and 2nd column for mask: dwi1,mask1
                                            dwi2,mask2 ...
        --ref_name VALUE:str                reference site name; required
        --resample VALUE:str                voxel size MxNxO to resample into
        --tar_list VALUE:ExistingFile       target csv/txt file with first column for dwi and 2nd column for mask: dwi1,mask1
                                            dwi2,mask2 ...
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
* [--nzero VALUE:str](#number-of-zero-padding)
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
You may test the whole pipeline as follows*:
    
    cd lib/tests
    ./test_pipeline.sh
    
NOTE: running the above test should take roughly an hour.

`./pipeline_test.sh` will download test data*, and run the whole processing pipeline on them. 
If the test is successful and complete, you should see the following output on the command line. 
    
    CONNECTOM mean FA:  0.675917931134666
    PRISMA mean FA before harmonization:  0.8008729791383536
    PRISMA mean FA after harmonization:  0.6366103024088171


\* If there is any problem downloading test data, try manually downloading and unzipping it to `lib/tests/` folder.

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

**NOTE** See [Caveats/Issues](#caveatsissues) that may occur while using many processors in parallel.

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


# Number of zero padding

During signal reconstruction, skull region can be noisy. So, each RISH feature is zero padded, median filtered in the 
skull region, and then used for signal reconstruction. The more noisy data you have, the more number of zeros you want 
to pad. The value should be 5-20.

    --nzero 10

# NRRD support

The pipeline is written for NIFTI image format. However, NRRD support is incorporated through [NIFTI --> NRRD](https://github.com/pnlbwh/Harmonization-Python/blob/parallel/lib/preprocess.py#L78) 
conversion on the fly.

See Tashrif Billah, Sylvain Bouix and Yogesh Rathi, Various MRI Conversion Tools, https://github.com/pnlbwh/conversion, 
2019, DOI: 10.5281/zenodo.2584003 for more details on the conversion method.


# Preprocessing

There are three preprocessing steps. More than one step can be specified at a time. However, they are performed in 
the following order.

## 1. Denoising
    
    --denoise        # turn on this flag to denoise voxel data

## 2. Bvalue mapping

    --bvalMap VALUE  # specify a bmax to scale bvalues into    

## 3. Resampling

    --resample VALUE # voxel size MxNxO to resample into


After preprocessing, the image lists are saved with `.modified` extension in the same location of provided lists, 
and used for further processing.
 
# Config   
    
With all the arguments passed to `harmonization.py` cli, a *config.ini* file is created in `lib/`. This file is used by 
submodules of the harmonization pipeline during executation. A sample `config.ini` file looks like the following:

    [DEFAULT]
    N_shm = 6
    N_proc = 8
    N_zero = 10
    resample = 1.5x1.5x1.5
    bvalMap = 1000
    denoise = 0
    travelHeads = 0
    debug = 1
    diffusionMeasures = MD,FA,GFA

# Template

Template directory is provided with `--template` argument. This directory is primarily used by `antsMultivariateTemplateConstruction2.sh` 
to save transform files. Then, various scales maps (*Scale_\*.nii.gz*) and mean templates (*Mean_\*.nii.gz*) 
are created in this directory. There are other template files: 
*Delta_\*.nii.gz*, *PercentageDiff_\*.nii.gz*, and *PercentageDiff_\*smooth.nii.gz*. You may look at them to know how good 
are the created templates.
    
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
    
RISH features and diffusion measures are computed from all the input images. Images of two modalites: 
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

During template construction, you want to provide a **small list of images** to `--ref_list` and `--tar_list`. 
[Section 2.4](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6230479/) of our paper gives a concise analysis about the number of samples to be chosen for template creation. The paper posits 20 as the "Gold Standard" for template creation while 16 as the minimum number. The rule of thumb should be to choose 16-20 subjects matched across sites for age, gender and IQ. However, 
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

* reference data is proprocessed, and registered to reference template space and then to MNI space ([IITmean_FA.nii.gz](https://www.nitrc.org/frs/download.php/6898/IITmean_FA.nii.gz))
* unprocessed target data is directly registered to MNI space
* harmonized target data is registered to target template space and then to MNI space
* once the data are in MNI space, we calculate mean FA over the [IITmean_FA_skeleton.nii.gz](https://www.nitrc.org/frs/download.php/6897/IITmean_FA_skeleton.nii.gz)
 
NOTE: Download the above data from [IIT HUMAN BRAIN ATLAS](http://www.iit.edu/~mri/IITHumanBrainAtlas.html) and place them in `IITAtlas/` directory.

The numbers should like like below:

    Printing statistics :
    CIDAR mean FA:  0.5217237675408243
    BSNIP mean FA before harmonization:  0.5072286796848892
    BSNIP mean FA after harmonization:  0.5321998242139347

As we see above, BSNIP (target) mean FA over [IITmean_FA_skeleton.nii.gz](https://www.nitrc.org/frs/download.php/6897/IITmean_FA_skeleton.nii.gz) 
for all the site images after harmonization increased to be almost equal to that of CIDAR (reference) mean FA.

Now there are two ways to debug:

## 1. With the pipeline 
Use `--debug` flag with any (or all) of `--create` and `--process`

## 2. Use seperately 
If you would like to debug at a later time, you need to specify three images lists:

* `--ref_list`: use the reference list with `.modified` extension
* `--tar_list`: use the unprocessed target list **without** the `.modified` extension
* `--harm_list`: use the harmonized target list that has `.harmonized` extension

    lib/harmonization.py --ref_list ref.txt.modified --tar_list target.csv --harm_list target.csv.modified.harmonized

NOTE: You should run the pipeline first before debugging separately because `--debug` makes use of files created 
in the pipeline.

## 3. With a list of FA images
The repository provides a more discrete discrete script for finding the goodness of harmonization. 

`$ lib/tests/fa_skeleton_test.py --help`
    
    usage: fa_skeleton_test.py [-h] -i INPUT -s SITE -t TEMPLATE
    
    Warps diffusion measures (FA, MD, GFA) to template space and then to subject
    space. Finally, calculates mean FA over IITmean_FA_skeleton.nii.gz
    
    optional arguments:
      -h, --help            show this help message and exit
      -i INPUT, --input INPUT
                            input list of FA images
      -s SITE, --site SITE  site name for locating template FA and mask in
                            tempalte directory
      -t TEMPLATE, --template TEMPLATE
                            template directory where Mean_{site}_FA.nii.gz and
                            {site}_Mask.nii.gz is located

This script does not depend of registration performed during the harmonization process. Rather, it performs all the 
steps mentioned above ([Debugging: In details](#debugging)) and computes mean FA over skeleton across all subjects 
in a site.

# Travel heads

Traveling heads mean *same set of subjects scanned with different scanners*. 
For example, scanning same set of people in Brigham and Women's Hospital and 
Massachusetts General Hospital. Therefore, users should set `--travelHeads` 
if they have same set of subjects across sites.


# Caveats/Issues

## 1. Template path

`antsMultivariateTemplateConstruction2.sh`: all the images need to have unique
prefix because transform files are created in the same `--template ./template/` directory. The value of `--template` 
should have `/` at the end. The pipeline appends one if there is not, but it is good to include it when specifying.

## 2. Multi-processing

[Multi threading](#-multi-threading) requires memory and processor availability. If pipeline does not continue past 
`unringing` or `shm_coeff` computation, your machine likely ran out of memory. Reducing `--nproc` to lower number of processors (i.e. 1-4) 
or swithcing to a powerful machine should help in this regard.

## 3. Tracker

In any case, feel free to submit an issue [here](https://github.com/pnlbwh/dMRIharmonization/issues). We shall get back to you as soon as possible.

# Reference

Zhang S, Arfanakis K. Evaluation of standardized and study-specific diffusion tensor imaging templates 
of the adult human brain: Template characteristics, spatial normalization accuracy, and detection of small 
inter-group FA differences. Neuroimage 2018;172:40-50.

Billah, Tashrif; Bouix, Sylvain; Rathi, Yogesh; Various MRI Conversion Tools, 
https://github.com/pnlbwh/conversion, 2019, DOI: 10.5281/zenodo.2584003.

