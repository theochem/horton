#!/usr/bin/env bash

source tools/qa/buildkite_common.sh
get_ancestor  # Writes $ANCESTOR_SHA variable.

echo "--- Basic source tests"
./tools/qa/check_names.py

echo "--- Build Cython files & HORTON"
./cleanfiles.sh
CC="ccache gcc" ./setup.py build_ext -i

echo "--- Packing build"
find horton -name "*.so" -o -name "*.pyc" | tar -zcvf horton_pr.tar.gz -T -
buildkite-agent artifact upload horton_pr.tar.gz

if [ "$BUILDKITE_PULL_REQUEST" != "false" ]; then
    echo "--- Checkout ancestor"
    git checkout $ANCESTOR_SHA

    echo "--- Build Cython files & HORTON [Ancestor]"
    ./cleanfiles.sh
    CC="ccache gcc" ./setup.py build_ext -i

    echo "--- Packing build [Ancestor]"
    find horton -name "*.so" -o -name "*.pyc" | tar -zcvf horton_ancestor.tar.gz -T -
    buildkite-agent artifact upload horton_ancestor.tar.gz
fi
