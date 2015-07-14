#!/usr/bin/env bash

# exit script when a line fails (instead of going to the next line)
set -e

COLOR='\033[0;32m'
RESET='\033[0m'

# check to see if h5py-2.2.1 is already present
if [ ! -d "$HOME/depinstall/h5py-2.2.1/lib/python2.7/site-packages" ]; then
(
    echo -e "${COLOR}Building and installing H5Py 2.2.1 from scratch${RESET}"
    cd
    mkdir -p depbuild
    cd depbuild
    wget http://pypi.python.org/packages/source/h/h5py/h5py-2.2.1.tar.gz
    tar -xzf h5py-2.2.1.tar.gz
    cd h5py-2.2.1
    echo "Actual build and install. This may take a while."
    python setup.py install --prefix=$HOME/depinstall/h5py-2.2.1 &> install.log
    tail install.log
)
else
    echo -e "${COLOR}Using Cached H5Py 2.2.1${RESET}"
fi
