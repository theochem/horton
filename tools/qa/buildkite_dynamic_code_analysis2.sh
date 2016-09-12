#!/usr/bin/env bash

if [ "$BUILDKITE_PULL_REQUEST" != "false" ]; then
    source tools/qa/buildkite_common.sh
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
    rm -rf $QAWORKDIR
    tools/qa/trapdoor_pylint.py feature

    echo "--- Unpack ancestor build from previous step"
    git checkout $ANCESTOR_SHA
    buildkite-agent artifact download horton_ancestor.tar.gz .
    ./cleanfiles.sh
    tar xvf horton_ancestor.tar.gz

    echo "--- Build refatoms"
    rm -rf data/refatoms/*.h5 data/refatoms/*.tar.bz2
    make -C data/refatoms/

    echo "--- Running trapdoor tests on ancestor branch"
    copy_qa_scripts

    $QAWORKDIR/trapdoor_pylint.py ancestor
    $QAWORKDIR/trapdoor_pylint.py report
fi
