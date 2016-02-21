#!/usr/bin/env bash

source tools/qa/common.sh

# This script assumes that a virtualenv is created, in which pip will then automatically
# install the dependencies.

# All packages are first wheeled, i.e. a binary package is made, which can be cached.

# Pip is called multiple times to have more regular output, which keeps travis happy.

pip install --upgrade pip
pip install wheel

# Packages for testing
pip wheel nose
pip install --upgrade nose
pip wheel pep8
pip install --upgrade pep8
pip wheel pylint
pip install --upgrade pylint
pip wheel coverage
pip install --upgrade coverage
pip wheel pep257
pip install --upgrade pep257

# Packages for HORTON
pip wheel numpy
pip install --upgrade numpy
pip wheel scipy
pip install --upgrade scipy
pip wheel matplotlib
pip install --upgrade matplotlib
pip wheel cython
pip install --upgrade cython
pip wheel h5py
pip install --upgrade h5py

# Packages for HORTON documentation
pip wheel sphinx
pip install --upgrade sphinx
pip wheel breathe
pip install --upgrade breathe

# Remove old wheels
./tools/qa/remove_old_wheels.py
