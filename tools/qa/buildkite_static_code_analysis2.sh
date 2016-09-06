#!/usr/bin/env bash

if [ "$BUILDKITE_PULL_REQUEST" != "false" ]; then
    source tools/qa/buildkite_common.sh
    checkout_merge_commit

    set +eo pipefail

    echo "--- Running trapdoors for feature branch"
    ./tools/qa/trapdoor_pylint.py feature || report_error "Trapdoor pylint failed (feature branch)"


    echo "--- Copying QA tools from feature branch"
    QAWORKDIR=`pwd`/qaworkdir
    cp -av ./tools/qa/trapdoor*.py $QAWORKDIR

    echo "--- Checking out ancestor branch"
    checkout_ancestor

    echo "--- Running trapdoors for ancestor branch"
    ${QAWORKDIR}/trapdoor_pylint.py ancestor || report_error "Trapdoor pylint failed (ancestor)"


    echo "--- Generating reports"
    ${QAWORKDIR}/trapdoor_pylint.py report || report_error "Trapdoor pylint regressions"


    if [ "$NUM_FAILED" -gt 0 ]; then
        echo -e "${RED}SOME TESTS FAILED (current branch)${RESET}"
        exit 1
    fi
    echo -e "${GREEN}ALL TESTS PASSED (current branch)${RESET}"
fi

exit 0
