#!/usr/bin/env bash

source tools/qa/common.sh
export NUM_FAILED=0

report_error() {
    echo -e "${RED}FAIL: ${1}${RESET}"
    ((NUM_FAILED++))
}

### 2) Testing in the ancestor with the master branch, unless we are currently on the master branch.

CURRENT_BRANCH=$(git rev-parse --abbrev-ref HEAD)
ANCESTOR_COMMIT=$(git merge-base master ${CURRENT_BRANCH})
CURRENT_COMMIT=$(git rev-parse HEAD)
if [ "${CURRENT_BRANCH}" != 'master' ] && [ "${CURRENT_COMMIT}" != ${ANCESTOR_COMMIT} ]; then
    # Copy the required scripts to the work directory, to make sure we're running with
    # the scripts from the feature branch, not the ancestor with the master branch.
    cp -av ./tools/qa/trapdoor*.py ${QAWORKDIR}/

    # Switch to ancestor
    git checkout ${ANCESTOR_COMMIT}

    # Clean stuff
    ./cleanfiles.sh &> /dev/null
    # Construct the reference atoms
    echo 'Rebuilding database of reference atoms'
    rm -rf data/refatoms/*.h5
    (cd data/refatoms; make all) || report_error "Failed to make reference atoms (ancestor)"
    # In-place build of HORTON
    python setup.py build_ext -i || report_error "Failed to build HORTON (ancestor)"

    # Run trapdoor tests (from QAWORKDIR)
    ${QAWORKDIR}/trapdoor_coverage.py ancestor || report_error "Trapdoor coverage failed (ancestor)"
    ${QAWORKDIR}/trapdoor_cppcheck.py ancestor || report_error "Trapdoor cppcheck failed (ancestor)"
    ${QAWORKDIR}/trapdoor_cpplint.py ancestor || report_error "Trapdoor cpplint failed (ancestor)"
    ${QAWORKDIR}/trapdoor_doxygen.py ancestor || report_error "Trapdoor doxygen failed (ancestor)"
    ${QAWORKDIR}/trapdoor_pylint.py ancestor || report_error "Trapdoor pylint failed (ancestor)"
    ${QAWORKDIR}/trapdoor_pycodestyle.py ancestor || report_error "Trapdoor pycodestyle failed (ancestor)"
    ${QAWORKDIR}/trapdoor_pydocstyle.py ancestor || report_error "Trapdoor pydocstyle failed (ancestor)"
    ${QAWORKDIR}/trapdoor_import.py ancestor || report_error "Trapdoor import failed (ancestor)"
    ${QAWORKDIR}/trapdoor_namespace.py ancestor || report_error "Trapdoor namespace failed (ancestor)"

    # Analyze trapdoor results (from QAWORKDIR)
    ${QAWORKDIR}/trapdoor_coverage.py report || report_error "Trapdoor coverage regressions"
    ${QAWORKDIR}/trapdoor_cppcheck.py report || report_error "Trapdoor cppcheck regressions"
    ${QAWORKDIR}/trapdoor_cpplint.py report || report_error "Trapdoor cpplint regressions"
    ${QAWORKDIR}/trapdoor_doxygen.py report || report_error "Trapdoor doxygen regressions"
    ${QAWORKDIR}/trapdoor_pylint.py report || report_error "Trapdoor pylint regressions"
    ${QAWORKDIR}/trapdoor_pycodestyle.py report || report_error "Trapdoor pycodestyle regressions"
    ${QAWORKDIR}/trapdoor_pydocstyle.py report || report_error "Trapdoor pydocstyle regressions"
    ${QAWORKDIR}/trapdoor_import.py report || report_error "Trapdoor import regressions"
    ${QAWORKDIR}/trapdoor_namespace.py report || report_error "Trapdoor namespace regressions"

    git checkout ${CURRENT_BRANCH}

    # Conclude
    if [ "$NUM_FAILED" -gt 0 ]; then
        echo -e "${RED}SOME TESTS FAILED, SEE ABOVE.${RESET}"
        exit 1
    fi

    echo -e "${GREEN}ALL TESTS PASSED (ancestor)${RESET}"
    exit 0
else
    echo "No need to test ancestor because HEAD is (merged in) master."
fi
