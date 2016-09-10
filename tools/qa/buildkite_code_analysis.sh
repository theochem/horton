#!/usr/bin/env bash

dir="tools/qa"
script="$dir/simulate_trapdoor_pr.py"

if [ "$BUILDKITE_PULL_REQUEST" != "false" ]; then
    source $dir/buildkite_common.sh
    get_ancestor  # Writes $ANCESTOR_SHA variable.

    echo "--- Prep working directory"
    rm -rf *_pr.tar.gz *_ancestor.tar.gz
    ./cleanfiles.sh

    echo "--- Build refatoms"
    rm -rf data/refatoms/*.h5 data/refatoms/*.tar.bz2
    make -C data/refatoms/

    echo "--- Running trapdoors tests"
    $dir/check_whitespace.py $ANCESTOR_SHA || report_error "Whitespace errors in some commits"

    export PYTHONPATH=`pwd`
    export HORTONDATA=`pwd`/data

    echo "--- Unpack PR build from previous step"
    buildkite-agent artifact download horton_pr.tar.gz .
    tar xvf horton_pr.tar.gz

    echo "--- Running trapdoor tests on PR branch"
    $script -vA $ANCESTOR_SHA -o feature $dir/trapdoor_coverage.py -t='--nproc=6'
    $script -vA $ANCESTOR_SHA -o feature $dir/trapdoor_namespace.py

    echo "--- Unpack ancestor build from previous step"
    buildkite-agent artifact download horton_ancestor.tar.gz .
    ./cleanfiles.sh
    tar xvf horton_ancestor.tar.gz

    echo "--- Build refatoms"
    rm -rf data/refatoms/*.h5 data/refatoms/*.tar.bz2
    make -C data/refatoms/

    echo "--- Running trapdoor tests on ancestor branch"
    for i in "ancestor" "report"; do
        $script -vA $ANCESTOR_SHA -o $i $dir/trapdoor_coverage.py -t='--nproc=6'
        $script -vA $ANCESTOR_SHA -o $i $dir/trapdoor_namespace.py
    done
fi
