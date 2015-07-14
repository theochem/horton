#!/usr/bin/env bash

# exit script when a line fails (instead of going to the next line)
set -e

COLOR='\033[0;32m'
RESET='\033[0m'

# check to see if numpy-1.9.1 is already present
if [ ! -d "$HOME/depinstall/numpy-1.9.1/lib/python2.7/site-packages" ]; then
(
    echo -e "${COLOR}Building and installing NumPy 1.9.1 from scratch${RESET}"
    cd
    mkdir -p depbuild
    cd depbuild
    wget http://pypi.python.org/packages/source/n/numpy/numpy-1.9.1.tar.gz
    tar -xzf numpy-1.9.1.tar.gz
    cd numpy-1.9.1
    echo "Actual build and install. This may take a while."
    python setup.py install --prefix=$HOME/depinstall/numpy-1.9.1 &> install.log
    tail install.log
)
else
    echo -e "${COLOR}Using Cached NumPy 1.9.1${RESET}"
fi
