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

echo "--- Running Nosetests (slow)"
nosetests -v --processes=2 --process-timeout=60 -A slow horton

echo "--- Running Nosetests (regression tests)"
nosetests -v --processes=2 --process-timeout=60 -A rt horton

if [ "$BUILDKITE_PULL_REQUEST" = "false" ]; then
  echo "--- Running Nosetests (fast)"
  nosetests -v --processes=2 --process-timeout=60 -A 'not (slow or rt)' horton
fi
