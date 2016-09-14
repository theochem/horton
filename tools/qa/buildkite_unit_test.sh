#!/usr/bin/env bash

source tools/qa/buildkite_common.sh

echo "--- Prep working directory"
rm -rf *_pr.tar.gz *_ancestor.tar.gz
./cleanfiles.sh

echo "--- Build refatoms"
rm -rf data/refatoms/*.h5 #data/refatoms/*.tar.bz2
make -C data/refatoms/

export PYTHONPATH=$PWD
export HORTONDATA=$PWD/data

echo "--- Unpack PR build from previous step"
buildkite-agent artifact download horton_pr.tar.gz .
tar xvf horton_pr.tar.gz

echo "--- Running Nosetests"
nosetests -v --processes=2 --process-timeout=60 -a slow horton

if [ "$BUILDKITE_PULL_REQUEST" = "false" ]; then
  nosetests -v --processes=2 --process-timeout=60 -a "!slow" horton
fi
