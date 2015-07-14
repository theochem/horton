#!/usr/bin/env bash

# exit script when a line fails (instead of going to the next line)
set -e

COLOR='\033[0;32m'
RESET='\033[0m'

# check to see if cython-0.17.1 is already present
if [ ! -d "$HOME/depinstall/cython-0.17.1/lib/python2.7/site-packages" ]; then
(
    echo -e "${COLOR}Building and installing Cython 0.17.1 from scratch${RESET}"
    cd
    mkdir -p depbuild
    cd depbuild
    wget https://pypi.python.org/packages/source/C/Cython/Cython-0.17.1.tar.gz
    tar -xzf Cython-0.17.1.tar.gz
    cd Cython-0.17.1
    echo "Actual build and install. This may take a while."
    python setup.py install --prefix=$HOME/depinstall/cython-0.17.1 &> install.log
    tail install.log
)
else
    echo -e "${COLOR}Using Cached Cython 0.17.1${RESET}"
fi
