#!/usr/bin/env bash

source tools/qa/common.sh
export NUM_FAILED=0

report_error() {
    echo -e "${RED}${1}${RESET}"
    ((NUM_FAILED++))
}

### Testing in the current branch

### a) Parts that are always done

# Check the author names
./tools/qa/check_names.py || report_error "Failed author/committer check (current branch)"
# Clean stuff
echo 'Cleaning source tree'
./cleanfiles.sh &> /dev/null
# Construct the reference atoms
echo 'Rebuilding database of reference atoms'
rm -rf data/refatoms/*.h5
(cd data/refatoms; make all) || report_error "Failed to make reference atoms (current branch)"
# In-place build of HORTON
python setup.py build_ext -i || report_error "Failed to build HORTON (current branch)"
# Run the slow tests
nosetests -v -a slow || report_error "Some slow tests failed (current branch)"
# Build the documentation
(cd doc; make html) || report_error "Failed to build documentation (current branch)"

### b) Parts that depend on the current branch
CURRENT_BRANCH=$(git rev-parse --abbrev-ref HEAD)
ANCESTOR_COMMIT=$(git merge-base master ${CURRENT_BRANCH})
CURRENT_COMMIT=$(git rev-parse HEAD)
if [ "${CURRENT_BRANCH}" == 'master' ] || [ "${CURRENT_COMMIT}" == ${ANCESTOR_COMMIT} ]; then
    # Run the fast tests
    nosetests -v -a '!slow' || report_error "Some fast tests failed (master branch)"
else
    # Check for whitespace errors in every commit.
    ./tools/qa/check_whitespace.py

    # Run the first part of the comparative tests.
    ./tools/qa/trapdoor_coverage.py feature || report_error "Trapdoor coverage failed (feature branch)"
    ./tools/qa/trapdoor_cppcheck.py feature || report_error "Trapdoor cppcheck failed (feature branch)"
    ./tools/qa/trapdoor_cpplint.py feature || report_error "Trapdoor cpplint failed (feature branch)"
    ./tools/qa/trapdoor_doxygen.py feature || report_error "Trapdoor doxygen failed (feature branch)"
    ./tools/qa/trapdoor_pylint.py feature || report_error "Trapdoor pylint failed (feature branch)"
    ./tools/qa/trapdoor_pycodestyle.py feature || report_error "Trapdoor pycodestyle failed (feature branch)"
    ./tools/qa/trapdoor_pydocstyle.py feature || report_error "Trapdoor pydocstyle failed (feature branch)"
    ./tools/qa/trapdoor_import.py feature || report_error "Trapdoor import failed (feature branch)"
    ./tools/qa/trapdoor_namespace.py feature || report_error "Trapdoor namespace failed (feature branch)"
fi

# Conclude
if [ "$NUM_FAILED" -gt 0 ]; then
    echo -e "${RED}SOME TESTS FAILED (current branch)${RESET}"
    exit 1
fi
echo -e "${GREEN}ALL TESTS PASSED (current branch)${RESET}"
exit 0
