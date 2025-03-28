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
    yum -y install wget file bzip2 vim git make unzip \
    libxcrypt-compat libXext libXt
    yum clean all
    rm -rf /var/cache/yum

    REPO=dMRIharmonization
    git clone --single-branch --branch shutil-and-bspline https://github.com/pnlbwh/$REPO.git

    # conda environment
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
    bash Miniconda3-latest-Linux-x86_64.sh -b -p miniconda3
    source miniconda3/bin/activate
    conda create -y -n harmonization -c conda-forge --override-channels python=3.8
    conda activate harmonization
    pip install -r $REPO/requirements.txt
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

    # fsl-6.0.7
    wget https://fsl.fmrib.ox.ac.uk/fsldownloads/fslinstaller.py -O fslinstaller.py > /dev/null 2>&1
    VERSION=6.0.7
    python fslinstaller.py -V $VERSION -d $HOME/fsl-$VERSION > /dev/null 2>&1
    FSLDIR=$HOME/fsl-$VERSION
    . $FSLDIR/etc/fslconf/fsl.sh
    
    # unring
    git clone https://bitbucket.org/reisert/unring.git

    # clean up sources
    conda deactivate
    conda clean -y --all
    rm -rf Miniconda3-latest-Linux-x86_64.sh .cache/pip/* $MCR.zip $MCR MCR_R2017a_Update_3_glnxa64.sh fslinstaller.py
    
    # provide write permissions
    chmod a+w $HOME
    

%environment
    # set up bashrc i.e shell
    export HOME=/home/pnlbwh/
    export USER=`whoami`
    export LANG=en_US.UTF-8

    FSLDIR=$HOME/fsl-6.0.7
    . ${FSLDIR}/etc/fslconf/fsl.sh
    
    MCRROOT=$HOME/MATLAB_Runtime/v92
    PATH=$HOME/unring/fsl/:${FSLDIR}/share/fsl/bin:$PATH
    CONDA_PREFIX=$HOME/miniconda3/envs/harmonization
    ANTSPATH=${CONDA_PREFIX}/bin

    PATH=${CONDA_PREFIX}/bin:$PATH

    export MCRROOT ANTSPATH PATH FSLDIR


