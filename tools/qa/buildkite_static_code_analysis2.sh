#!/usr/bin/env bash

dir="tools/qa"

if [ "$BUILDKITE_PULL_REQUEST" != "false" ]; then
    source $dir/buildkite_common.sh
    get_ancestor  # Writes $ANCESTOR_SHA variable.

    echo "--- Prep working directory"
    rm -rf horton_pr.tar.gz horton_ancestor.tar.gz
    ./cleanfiles.sh

    echo "--- Unpack PR build from previous step"
    buildkite-agent artifact download horton_pr.tar.gz .
    tar xvf horton_pr.tar.gz

    echo "--- Running trapdoors tests on PRbranch"
    $dir/simulate_trapdoor_pr.py -vA $ANCESTOR_SHA -o feature $dir/trapdoor_pylint.py

    echo "--- Unpack ancestor build from previous step"
    buildkite-agent artifact download horton_ancestor.tar.gz .
    ./cleanfiles.sh
    tar xvf horton_ancestor.tar.gz

    echo "--- Running trapdoor tests on ancestor branch"
    for i in "ancestor" "report"; do
        $dir/simulate_trapdoor_pr.py -vA $ANCESTOR_SHA -o $i $dir/trapdoor_pylint.py
    done
fi
