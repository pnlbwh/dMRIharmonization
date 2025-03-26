Bootstrap: docker
From: redhat/ubi9:9.5-1738643550

%labels
    MAINTAINER Tashrif Billah <tbillah@bwh.harvard.edu>

%help
    For running container
        - https://github.com/pnlbwh/dMRIharmonization

    Please report issues on GitHub.


%post
    # set up user and working directory
    mkdir /home/pnlbwh
    cd /home/pnlbwh
    export HOME=`pwd`

    # install required libraries
    yum -y install wget file bzip2 vim git make unzip
    yum clean all
    rm -rf /var/cache/yum

    git clone --single-branch --branch shutil-and-bspline https://github.com/pnlbwh/dMRIharmonization.git

    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
    bash Miniconda3-latest-Linux-x86_64.sh -b -p miniconda3
    source miniconda3/bin/activate
    conda create -y -n harmonization -c conda-forge --override-channels python=3.8
    conda activate harmonization
    pip install -r dMRIharmonization/requirements.txt
    conda install -y pnlbwh::ants

    # MCR 2017a
    MCR=MCR_R2017a_glnxa64_installer
    wget https://ssd.mathworks.com/supportfiles/downloads/R2017a/deployment_files/R2017a/installers/glnxa64/$MCR.zip
    unzip $MCR.zip -d $MCR
    cd $MCR
    ./install -mode silent -agreeToLicense yes -destinationFolder $HOME/MATLAB_Runtime
    cd ..

    # MCR 2017a security updates
    wget https://ssd.mathworks.com/supportfiles/downloads/R2017a/deployment_files/R2017a/installers/glnxa64/MCR_R2017a_Update_3_glnxa64.sh
    bash MCR_R2017a_Update_3_glnxa64.sh -d=`pwd`/MATLAB_Runtime/v92
    

%environment
    # set up bashrc i.e shell
    export HOME=/home/pnlbwh/
    export USER=`whoami`
    export LANG=en_US.UTF-8

    NEW_SOFT_DIR=$HOME
    source ${NEW_SOFT_DIR}/fsl-6.0.7/activate.sh
    MCRROOT=${NEW_SOFT_DIR}/MATLAB_Runtime/v92
    PATH=${NEW_SOFT_DIR}/unring/fsl/:$PATH
    CONDA_PREFIX=${NEW_SOFT_DIR}/miniconda3/envs/harmonization
    ANTSPATH=${CONDA_PREFIX}/bin

    PATH=${CONDA_PREFIX}/bin:$PATH

    export MCRROOT ANTSPATH PATH


