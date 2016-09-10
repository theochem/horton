#!/usr/bin/env bash

dir="tools/qa"
script="$dir/simulate_trapdoor_pr.py"

if [ "$BUILDKITE_PULL_REQUEST" != "false" ]; then
    source $dir/buildkite_common.sh
    get_ancestor  # Writes $ANCESTOR_SHA variable.

    echo "--- Prep working directory"
    rm -rf *_pr.tar.gz *_ancestor.tar.gz
    ./cleanfiles.sh

    echo "--- Load PR refatoms"
    buildkite-agent artifact download refatoms_pr.tar.gz .
    tar xvf refatoms_pr.tar.gz

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

    echo "--- Load ancestor refatoms"
    buildkite-agent artifact download refatoms_ancestor.tar.gz .
    rm -rf data/refatoms/*
    tar xvf refatoms_ancestor.tar.gz

    echo "--- Running trapdoor tests on ancestor branch"
    for i in "ancestor" "report"; do
        $script -vA $ANCESTOR_SHA -o $i $dir/trapdoor_coverage.py -t='--nproc=6'
        $script -vA $ANCESTOR_SHA -o $i $dir/trapdoor_namespace.py
    done
fi
