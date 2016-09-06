#!/usr/bin/env bash

source tools/qa/buildkite_common.sh
checkout_merge_commit

echo "--- Downloading build artifact"
#buildkite-agent artifact download horton-install.tar.gz
tar -xvf horton-install.tar.gz

echo "--- Running Nosetests"
cd installation
PATH=$PATH:`pwd`/bin PYTHONPATH=`pwd`/lib/python2.7/site-packages:`pwd`/lib64/python2.7/site-packages HORTONDATA=`pwd`/share/horton nosetests -v --processes=2 --process-timeout=60 -a slow horton

if [ "$BUILDKITE_PULL_REQUEST" = "false" ]; then
  PATH=$PATH:`pwd`/bin PYTHONPATH=`pwd`/lib/python2.7/site-packages:`pwd`/lib64/python2.7/site-packages HORTONDATA=`pwd`/share/horton nosetests -v --processes=2 --process-timeout=60 -a "!slow" horton
fi