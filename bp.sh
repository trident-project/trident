# This shell script installs both Trident as well as its dependencies.
# It is needed by Bitbucket's Pipelines functionality for continuous 
# integration of automated testing.

# install conda and yt dependencies with conda
cd ..
wget --quiet https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh
bash ./Miniconda2-latest-Linux-x86_64.sh -b -p ./trident-conda -f
export PATH=$PWD/trident-conda/bin:$PATH
conda install -q -y mercurial cython h5py matplotlib sympy numpy pytest flake8

# set up development build of yt
hg clone https://bitbucket.org/yt_analysis/yt yt-hg
cd yt-hg
hg up yt
pip install -e .
cd ..

# repo has been cloned here
cd $BITBUCKET_CLONE_DIR
hg up tip
pip install -e .

# install data file and create config.tri file in CWD
# so trident works, despite not having access to ~/.trident
cd tests
export RUN_DOWNLOAD_TEST=1
wget http://trident-project.org/data/ion_table/config.tri
wget http://trident-project.org/data/ion_table/hm2012_lr.h5.gz
gunzip hm2012_lr.h5.gz
py.test test_download.py
# start the tests themselves
export RUN_DOWNLOAD_TEST=0
py.test
