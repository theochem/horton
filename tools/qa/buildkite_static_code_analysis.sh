#!/usr/bin/env bash

if [ "$BUILDKITE_PULL_REQUEST" != "false" ]; then
    source tools/qa/buildkite_common.sh
    get_ancestor  # Writes $ANCESTOR_SHA variable.

    echo "--- Prep working directory"
    ./cleanfiles.sh

    echo "--- Checking for whitespace errors in every commit"
    tools/qa/check_commits.py $ANCESTOR_SHA || report_error "Bad habits detected in some commits"

    PATH=$PATH:~/.local/bin  # fix for ubuntu paths

    echo "--- Running trapdoors tests on PR branch"
    rm -rf $QAWORKDIR/*.pp

    TRAPDOORS="trapdoor_cppcheck.py
    trapdoor_cpplint.py
    trapdoor_doxygen.py
    trapdoor_import.py
    trapdoor_pycodestyle.py
    trapdoor_pydocstyle.py"

    for i in ${TRAPDOORS}; do
        tools/qa/$i feature
    done

    echo "--- Copy PR version of trapdoor scripts to QAWORKDIR"
    copy_trapdoor_scripts

    echo "--- Checkout ancestor, rerun trapdoor scripts and report"
    git checkout $ANCESTOR_SHA

    for i in ${TRAPDOORS}; do
        $QAWORKDIR/$i ancestor
    done

    for i in ${TRAPDOORS}; do
        $QAWORKDIR/$i report
    done
fi
