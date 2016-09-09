#!/usr/bin/env bash

source tools/qa/buildkite_common.sh
checkout_merge_commit

echo "--- Basic source tests"
./tools/qa/check_names.py

echo "--- Build refatoms"
rm -rf data/refatoms/*.h5 data/refatoms/*.tar.bz2
make -C data/refatoms/

echo "--- Build Cython files & HORTON"
./cleanfiles.sh
rm -rf installation
./setup.py install --prefix=`pwd`/installation

echo "--- Running Nosetests"
# Configure nose to run in parallel. Not sure if this is a good idea.
printf "[nosetests]\nprocesses=4\nprocess-timeout=600\n" > ~/.noserc


cd installation
#PATH=$PATH:`pwd`/bin PYTHONPATH=`pwd`/lib/python2.7/site-packages:`pwd`/lib64/python2.7/site-packages HORTONDATA=`pwd`/share/horton nosetests -v -a slow horton
PATH=$PATH:`pwd`/bin PYTHONPATH=`pwd`/lib/python2.7/site-packages:`pwd`/lib64/python2.7/site-packages HORTONDATA=`pwd`/share/horton nosetests -v -a slow horton

if [ "$BUILDKITE_PULL_REQUEST" = "false" ]; then
  PATH=$PATH:`pwd`/bin PYTHONPATH=`pwd`/lib/python2.7/site-packages:`pwd`/lib64/python2.7/site-packages HORTONDATA=`pwd`/share/horton nosetests -v -a "!slow" horton
fi
