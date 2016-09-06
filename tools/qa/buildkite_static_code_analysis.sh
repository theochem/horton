#!/usr/bin/env bash

if [ "$BUILDKITE_PULL_REQUEST" != "false" ]; then
    source tools/qa/buildkite_common.sh
    checkout_merge_commit

    set +eo pipefail

    ORIG_PATH=$PATH

    export PATH=$ORIG_PATH:`pwd`/installation/bin
    export PYTHONPATH=`pwd`/installation/lib/python2.7/site-packages
    export PYTHONPATH=$PYTHONPATH:`pwd`/installation/lib64/python2.7/site-packages
    export HORTONDATA=`pwd`/installation/share/horton

    echo "--- Running coverage for feature branch"
    ./tools/qa/trapdoor_cppcheck.py feature || report_error "Trapdoor cppcheck failed (feature branch)"
    ./tools/qa/trapdoor_cpplint.py feature || report_error "Trapdoor cpplint failed (feature branch)"
    ./tools/qa/trapdoor_doxygen.py feature || report_error "Trapdoor doxygen failed (feature branch)"
    ./tools/qa/trapdoor_import.py feature || report_error "Trapdoor import failed (feature branch)"



    echo "--- Copying QA tools from feature branch"
    QAWORKDIR=`pwd`/qaworkdir
    cp -av ./tools/qa/trapdoor*.py $QAWORKDIR

    echo "--- Checking out ancestor branch"
    checkout_ancestor

    echo "--- Downloading build artifacts"
    buildkite-agent artifact download horton-ancestor-install.tar.gz
    tar -xvf horton-ancestor-install.tar.gz

    echo "--- Running coverage for ancestor branch"
    export PATH=$ORIG_PATH:`pwd`/ancestor_installation/bin
    export PYTHONPATH=`pwd`/ancestor_installation/lib/python2.7/site-packages
    export PYTHONPATH=$PYTHONPATH:`pwd`/ancestor_installation/lib64/python2.7/site-packages
    export HORTONDATA=`pwd`/ancestor_installation/share/horton

    ${QAWORKDIR}/trapdoor_cppcheck.py ancestor || report_error "Trapdoor cppcheck failed (ancestor)"
    ${QAWORKDIR}/trapdoor_cpplint.py ancestor || report_error "Trapdoor cpplint failed (ancestor)"
    ${QAWORKDIR}/trapdoor_doxygen.py ancestor || report_error "Trapdoor doxygen failed (ancestor)"
    ${QAWORKDIR}/trapdoor_import.py ancestor || report_error "Trapdoor import failed (ancestor)"


    echo "--- Generating reports"
    ${QAWORKDIR}/trapdoor_cppcheck.py report || report_error "Trapdoor cppcheck regressions"
    ${QAWORKDIR}/trapdoor_cpplint.py report || report_error "Trapdoor cpplint regressions"
    ${QAWORKDIR}/trapdoor_doxygen.py report || report_error "Trapdoor doxygen regressions"
    ${QAWORKDIR}/trapdoor_import.py report || report_error "Trapdoor import regressions"



    if [ "$NUM_FAILED" -gt 0 ]; then
        echo -e "${RED}SOME TESTS FAILED (current branch)${RESET}"
        exit 1
    fi
    echo -e "${GREEN}ALL TESTS PASSED (current branch)${RESET}"
fi

exit 0