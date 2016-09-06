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

    echo "--- Running trapdoors for feature branch"
    ./tools/qa/check_whitespace.py || report_error "Whitespace errors in some commits"
    ./tools/qa/trapdoor_coverage.py feature || report_error "Trapdoor coverage failed (feature branch)"
    ./tools/qa/trapdoor_namespace.py feature || report_error "Trapdoor namespace failed (feature branch)"

    echo "--- Copying QA tools from feature branch"
    QAWORKDIR=`pwd`/qaworkdir
    cp -av ./tools/qa/trapdoor*.py $QAWORKDIR

    echo "--- Checking out ancestor branch"
    checkout_ancestor

    echo "--- Downloading build artifacts"
    buildkite-agent artifact download horton-ancestor-install.tar.gz
    tar -xvf horton-ancestor-install.tar.gz

    echo "--- Running trapdoors for ancestor branch"
    export PATH=$ORIG_PATH:`pwd`/ancestor_installation/bin
    export PYTHONPATH=`pwd`/ancestor_installation/lib/python2.7/site-packages
    export PYTHONPATH=$PYTHONPATH:`pwd`/ancestor_installation/lib64/python2.7/site-packages
    export HORTONDATA=`pwd`/ancestor_installation/share/horton

    $QAWORKDIR/trapdoor_coverage.py ancestor || report_error "Trapdoor coverage failed (ancestor)"
    ${QAWORKDIR}/trapdoor_namespace.py ancestor || report_error "Trapdoor namespace failed (ancestor)"

    echo "--- Generating reports"
    $QAWORKDIR/trapdoor_coverage.py report || report_error "Trapdoor coverage regressions"
    ${QAWORKDIR}/trapdoor_namespace.py report || report_error "Trapdoor namespace regressions"

    if [ "$NUM_FAILED" -gt 0 ]; then
        echo -e "${RED}SOME TESTS FAILED (current branch)${RESET}"
        exit 1
    fi
    echo -e "${GREEN}ALL TESTS PASSED (current branch)${RESET}"
fi

exit 0