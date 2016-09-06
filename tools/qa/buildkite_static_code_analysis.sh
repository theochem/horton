#!/usr/bin/env bash

if [ "$BUILDKITE_PULL_REQUEST" != "false" ]; then
    source tools/qa/buildkite_common.sh
    checkout_merge_commit

    set +eo pipefail

    echo "--- Running trapdoors for feature branch"
    ./tools/qa/trapdoor_cppcheck.py feature || report_error "Trapdoor cppcheck failed (feature branch)"
    ./tools/qa/trapdoor_cpplint.py feature || report_error "Trapdoor cpplint failed (feature branch)"
    ./tools/qa/trapdoor_doxygen.py feature || report_error "Trapdoor doxygen failed (feature branch)"
    ./tools/qa/trapdoor_import.py feature || report_error "Trapdoor import failed (feature branch)"
    ./tools/qa/trapdoor_pycodestyle.py feature || report_error "Trapdoor pycodestyle failed (feature branch)"
    ./tools/qa/trapdoor_pydocstyle.py feature || report_error "Trapdoor pydocstyle failed (feature branch)"



    echo "--- Copying QA tools from feature branch"
    QAWORKDIR=`pwd`/qaworkdir
    cp -av ./tools/qa/trapdoor*.py $QAWORKDIR

    echo "--- Checking out ancestor branch"
    checkout_ancestor

    echo "--- Running trapdoors for ancestor branch"
    ${QAWORKDIR}/trapdoor_cppcheck.py ancestor || report_error "Trapdoor cppcheck failed (ancestor)"
    ${QAWORKDIR}/trapdoor_cpplint.py ancestor || report_error "Trapdoor cpplint failed (ancestor)"
    ${QAWORKDIR}/trapdoor_doxygen.py ancestor || report_error "Trapdoor doxygen failed (ancestor)"
    ${QAWORKDIR}/trapdoor_import.py ancestor || report_error "Trapdoor import failed (ancestor)"
    ${QAWORKDIR}/trapdoor_pycodestyle.py ancestor || report_error "Trapdoor pycodestyle failed (ancestor)"
    ${QAWORKDIR}/trapdoor_pydocstyle.py ancestor || report_error "Trapdoor pydocstyle failed (ancestor)"


    echo "--- Generating reports"
    ${QAWORKDIR}/trapdoor_cppcheck.py report || report_error "Trapdoor cppcheck regressions"
    ${QAWORKDIR}/trapdoor_cpplint.py report || report_error "Trapdoor cpplint regressions"
    ${QAWORKDIR}/trapdoor_doxygen.py report || report_error "Trapdoor doxygen regressions"
    ${QAWORKDIR}/trapdoor_import.py report || report_error "Trapdoor import regressions"
    ${QAWORKDIR}/trapdoor_pycodestyle.py report || report_error "Trapdoor pycodestyle regressions"
    ${QAWORKDIR}/trapdoor_pydocstyle.py report || report_error "Trapdoor pydocstyle regressions"



    if [ "$NUM_FAILED" -gt 0 ]; then
        echo -e "${RED}SOME TESTS FAILED (current branch)${RESET}"
        exit 1
    fi
    echo -e "${GREEN}ALL TESTS PASSED (current branch)${RESET}"
fi

exit 0