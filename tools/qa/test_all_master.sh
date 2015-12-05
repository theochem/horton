#!/usr/bin/env bash

source tools/qa/deps/common.sh
export NUM_FAILED=0

report_error() {
    echo -e "${RED}FAIL: ${1}${RESET}"
    ((NUM_FAILED++))
}

abort_error() {
    echo -e "${RED}${1}${RESET}"
    echo -e "${RED}TESTS ABORTED (master branch)${RESET}"
    exit 1
}


### 2) Testing in the master branch (if that isn't current)

CURRENT_BRANCH=$(git rev-parse --abbrev-ref HEAD)
if [ "${CURRENT_BRANCH}" != 'master' ]; then
    # Check if the feature branch is a proper descendant of the master branch. If not,
    # abort the tests.
    git merge-base --is-ancestor master ${CURRENT_BRANCH} || abort_error "Feature branch is not a direct descendant of master."

    # Copy the required scripts to the work directory, to make sure we're running with
    # the scripts from the feature branch, not the master branch
    cp -av ./tools/qa/trapdoor*.py ${QAWORKDIR}/

    # Switch to master
    git checkout master

    # Needed for coverage: rebuild
    # Activate dependencies (from master branch, may be different)
    for DEPDIR in $(cat tools/qa/deps/dirs.txt); do
        [[ -f "tools/qa/deps/${DEPDIR}/activate.sh" ]] && source tools/qa/deps/${DEPDIR}/activate.sh
    done
    # Clean stuff
    ./cleanfiles.sh &> /dev/null
    rm -rf data/refatoms/*.h5
    # Construct the reference atoms
    (cd data/refatoms; make all) || report_error "Failed to make reference atoms (master branch)"
    # In-place build of HORTON
    python setup.py build_ext -i -L ${LD_LIBRARY_PATH} || report_error "Failed to build HORTON (master branch)"

    # Run trapdoor tests (from QAWORKDIR)
    ${QAWORKDIR}/trapdoor_coverage.py master || report_error "Trapdoor coverage failed (master branch)"
    ${QAWORKDIR}/trapdoor_cppcheck.py master || report_error "Trapdoor cppcheck failed (master branch)"
    ${QAWORKDIR}/trapdoor_cpplint.py master || report_error "Trapdoor cpplint failed (master branch)"
    ${QAWORKDIR}/trapdoor_doxygen.py master || report_error "Trapdoor doxygen failed (master branch)"
    ${QAWORKDIR}/trapdoor_pylint.py master || report_error "Trapdoor pylint failed (master branch)"
    ${QAWORKDIR}/trapdoor_pep8.py master || report_error "Trapdoor pep8 failed (master branch)"

    # Analyze trapdoor results (from QAWORKDIR)
    ${QAWORKDIR}/trapdoor_coverage.py report || report_error "Trapdoor coverage regressions"
    ${QAWORKDIR}/trapdoor_cppcheck.py report || report_error "Trapdoor cppcheck regressions"
    ${QAWORKDIR}/trapdoor_cpplint.py report || report_error "Trapdoor cpplint regressions"
    ${QAWORKDIR}/trapdoor_doxygen.py report || report_error "Trapdoor doxygen regressions"
    ${QAWORKDIR}/trapdoor_pylint.py report || report_error "Trapdoor pylint regressions"
    ${QAWORKDIR}/trapdoor_pep8.py report || report_error "Trapdoor pep8 regressions"

    git checkout ${CURRENT_BRANCH}

    # Conclude
    if [ "$NUM_FAILED" -gt 0 ]; then
        echo -e "${RED}SOME TESTS FAILED, SEE ABOVE.${RESET}"
        exit 1
    fi
fi

echo -e "${GREEN}ALL TESTS PASSED (master branch)${RESET}"
exit 0
