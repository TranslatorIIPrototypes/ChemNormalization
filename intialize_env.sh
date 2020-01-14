#!/bin/sh

# make a root directory for code
mkdir ./code

# go to the directory
cd code

# clone the cod directory
git clone https://github.com/TranslatorIIPrototypes/ChemNormalization.git

# go to the Chem Normalization directory
cd ChemNormalization

# get the anaconda install script
wget http://repo.continuum.io/archive/Anaconda3-4.0.0-Linux-x86_64.sh

# run the install script
bash Anaconda3-4.0.0-Linux-x86_64.sh

# update the now installed conda to the latest version
conda update

# add the conda-forge package install channel
conda config --add channels conda-forge

# Create a conda virtual environment
conda create -n chemNodeVenv python=3.7

# Install package requirements
conda install --yes --file requirements.txt
