#!/usr/bin/env bash

# exit script when a line fails (instead of going to the next line)
set -e

COLOR='\033[0;32m'
RESET='\033[0m'

# check to see if sympy-0.7.6 is already present
if [ ! -d "$HOME/depinstall/sympy-0.7.6/lib/python2.7/site-packages" ]; then
(
    echo -e "${COLOR}Building and installing SymPy 0.7.6 from scratch${RESET}"
    cd
    mkdir -p depbuild
    cd depbuild
    wget http://pypi.python.org/packages/source/s/sympy/sympy-0.7.6.tar.gz
    tar -xzf sympy-0.7.6.tar.gz
    cd sympy-0.7.6
    echo "Actual build and install. This may take a while."
    python setup.py install --prefix=$HOME/depinstall/sympy-0.7.6 &> install.log
    tail install.log
)
else
    echo -e "${COLOR}Using Cached SymPy 0.7.6${RESET}"
fi
