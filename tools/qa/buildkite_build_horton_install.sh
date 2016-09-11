#!/usr/bin/env bash

source tools/qa/buildkite_common.sh

echo "--- Basic source tests"
./tools/qa/check_names.py

echo "--- Build refatoms"
rm -rf data/refatoms/*.h5 data/refatoms/*.tar.bz2
make -C data/refatoms/

echo "--- Build Cython files & HORTON"
./cleanfiles.sh
rm -rf installation
./setup.py install --prefix=$PWD/installation

echo "--- Running Nosetests"
cd installation
#PATH=$PATH:$PWD/bin PYTHONPATH=$PWD/lib/python2.7/site-packages:$PWD/lib64/python2.7/site-packages HORTONDATA=$PWD/share/horton nosetests -v -a slow horton
PATH=$PATH:$PWD/bin PYTHONPATH=$PWD/lib/python2.7/site-packages:$PWD/lib64/python2.7/site-packages HORTONDATA=$PWD/share/horton nosetests -v --processes=2 --process-timeout=60 -a slow horton

if [ "$BUILDKITE_PULL_REQUEST" = "false" ]; then
  PATH=$PATH:$PWD/bin PYTHONPATH=$PWD/lib/python2.7/site-packages:$PWD/lib64/python2.7/site-packages HORTONDATA=$PWD/share/horton nosetests -v --processes=2 --process-timeout=60 -a "!slow" horton
fi
