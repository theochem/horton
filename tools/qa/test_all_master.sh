#!/usr/bin/env bash

source tools/qa/deps/common.sh
export NUM_FAILED=0

report_error() {
    echo -e "${RED}FAIL: ${1}${RESET}"
    ((NUM_FAILED++))
}


### 2) Testing in the master branch (if that isn't current)

CURRENT_BRANCH=$(git rev-parse --abbrev-ref HEAD)
if [ "${CURRENT_BRANCH}" != 'master' ]; then
    git checkout master

    # Needed for coverage: rebuild
    # Activate dependencies
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

    # Run trapdoor tests
    ./tools/qa/trapdoor_coverage.py master || report_error "Trapdoor coverage failed (master branch)"
    ./tools/qa/trapdoor_cppcheck.py master || report_error "Trapdoor cppcheck failed (master branch)"
    ./tools/qa/trapdoor_pylint.py master || report_error "Trapdoor pylint failed (master branch)"
    ./tools/qa/trapdoor_pep8.py master || report_error "Trapdoor pep8 failed (master branch)"

    # Analyze trapdoor results
    ./tools/qa/trapdoor_coverage.py report || report_error "Trapdoor coverage regressions (master branch)"
    ./tools/qa/trapdoor_cppcheck.py report || report_error "Trapdoor cppcheck regressions (master branch)"
    ./tools/qa/trapdoor_pylint.py report || report_error "Trapdoor pylint regressions (master branch)"
    ./tools/qa/trapdoor_pep8.py report || report_error "Trapdoor pep8 regressions (master branch)"

    git checkout ${CURRENT_BRANCH}

    # Conclude
    if [ "$NUM_FAILED" -gt 0 ]; then
        echo -e "${RED}SOME TESTS FAILED, SEE ABOVE.${RESET}"
        exit 1
    fi
fi

echo -e "${GREEN}ALL TESTS PASSED (master branch)${RESET}"
exit 0
