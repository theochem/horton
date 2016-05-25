#!/usr/bin/env bash

# This script does a complete local run-through of what runs on travis-ci.
# The following must be installed prior to running this script:
#   - virtualenv
#   - gcc
#   - g++
#   - gfortran
#   - libatlas-base-dev
#   - libhdf5-serial-dev
#   - libtool
#   - doxygen
#   - cppcheck

source tools/qa/common.sh
set -e

# 1) make a virtual environment and start using it
virtualenv ${QAWORKDIR}/venv
source ${QAWORKDIR}/venv/bin/activate

# 2) Set some environment variables for pip
export PIP_WHEEL_DIR=${CACHED}/wheels
export PIP_FIND_LINKS=file://${CACHED}/wheels
export PIP_NO_CACHE_DIR=false

# 3) Install dependencies in the virtualenv, just as on travis-ci
./tools/qa/install_deps_pip.sh

# 4) Install extra dependencies
./tools/qa/install_deps_extra_twobranches.sh

# 5) Run the actual tests (includes building HORTON).
./tools/qa/test_all_twobranches.sh
