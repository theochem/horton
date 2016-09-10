#!/usr/bin/env bash

dir="tools/qa"

if [ "$BUILDKITE_PULL_REQUEST" != "false" ]; then
    source $dir/buildkite_common.sh
    get_ancestor  # Writes $ANCESTOR_SHA variable.

    echo "--- Prep working directory"
    rm -rf *_pr.tar.gz *_ancestor.tar.gz
    ./cleanfiles.sh

    echo "--- Unpack PR build from previous step"
    buildkite-agent artifact download horton_pr.tar.gz .
    tar xvf horton_pr.tar.gz

    echo "--- Build refatoms"
    rm -rf data/refatoms/*.h5 data/refatoms/*.tar.bz2
    make -C data/refatoms/

    echo "--- Running trapdoors tests on PR branch"
    $dir/simulate_trapdoor_pr.py -vA $ANCESTOR_SHA -o feature $dir/trapdoor_pylint.py

    echo "--- Unpack ancestor build from previous step"
    buildkite-agent artifact download horton_ancestor.tar.gz .
    ./cleanfiles.sh
    tar xvf horton_ancestor.tar.gz

    echo "--- Build refatoms"
    rm -rf data/refatoms/*.h5 data/refatoms/*.tar.bz2
    make -C data/refatoms/

    echo "--- Running trapdoor tests on ancestor branch"
    for i in "ancestor" "report"; do
        $dir/simulate_trapdoor_pr.py -vA $ANCESTOR_SHA -o $i $dir/trapdoor_pylint.py
    done
fi
