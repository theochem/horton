#!/usr/bin/env bash

# exit script when a line fails (instead of going to the next line)
set -e

COLOR='\033[0;32m'
RESET='\033[0m'

# check to see if scipy-0.11.0 is already present
if [ ! -d "$HOME/depinstall/scipy-0.11.0/lib/python2.7/site-packages" ]; then
(
    echo -e "${COLOR}Building and installing SciPy 0.11.0 from scratch${RESET}"
    cd
    mkdir -p depbuild
    cd depbuild
    wget http://pypi.python.org/packages/source/s/scipy/scipy-0.11.0.tar.gz
    tar -xzf scipy-0.11.0.tar.gz
    cd scipy-0.11.0
    echo "Actual build and install. This may take a while."
    python setup.py install --prefix=$HOME/depinstall/scipy-0.11.0 &> install.log
    tail install.log
)
else
    echo -e "${COLOR}Using Cached SciPy 0.11.0${RESET}"
fi
