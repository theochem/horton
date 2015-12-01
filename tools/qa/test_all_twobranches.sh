#!/usr/bin/env bash

source tools/qa/deps/common.sh

### 1) Testing in the current branch
# This opens a new subprocess in which the dependencies of the current branch are activated.
(

### 1a) Parts that are always done
# Activate dependencies
for DEPDIR in $(cat tools/qa/deps/dirs.txt); do
    [[ -f "tools/qa/deps/${DEPDIR}/activate.sh" ]] && source tools/qa/deps/${DEPDIR}/activate.sh
done
## Clean stuff
#./cleanfiles.sh &> /dev/null
#rm -rf data/refatoms/*.h5
## Construct the reference atoms
#(cd data/refatoms; make all)
## In-place build of HORTON
#python setup.py build_ext -i -L ${LD_LIBRARY_PATH}
## Run the fast tests
#nosetests -v -a '!slow'
## Run the slow tests
#nosetests -v -a slow
## Build the documentation
#(cd doc; make html)

### 1b) Parts that depend on the current branch
CURRENT_BRANCH=$(git rev-parse --abbrev-ref HEAD)
if [ "${CURRENT_BRANCH}" != 'master' ]; then
    # Run the first part of the comparative tests if the current branch is not the
    # master branch.
    ./tools/qa/trapdoor_cppcheck.py feature
    ./tools/qa/trapdoor_pylint.py feature
    ./tools/qa/trapdoor_pep8.py feature
fi

)


### 2) Testing in the master branch (if that isn't current)
# This opens a new subprocess in which we will activate the dependencies of the current
# branch.
(

CURRENT_BRANCH=$(git rev-parse --abbrev-ref HEAD)
if [ "${CURRENT_BRANCH}" != 'master' ]; then
    git checkout master

    # Needed for coverage: rebuild
    ## Activate dependencies
    #for DEPDIR in $(cat tools/qa/deps/dirs.txt); do
    #    [[ -f "tools/qa/deps/${DEPDIR}/activate.sh" ]] && source tools/qa/deps/${DEPDIR}/activate.sh
    #done
    #done
    ## Clean stuff
    #./cleanfiles.sh &> /dev/null
    #rm -rf data/refatoms/*.h5
    ## Construct the reference atoms
    #(cd data/refatoms; make all)
    ## In-place build of HORTON
    #python setup.py build_ext -i -L ${LD_LIBRARY_PATH}

    ./tools/qa/trapdoor_cppcheck.py master
    ./tools/qa/trapdoor_pylint.py master
    ./tools/qa/trapdoor_pep8.py master

    NUM_FAILED=0
    ./tools/qa/trapdoor_cppcheck.py report
    [[ "$?" != 0 ]] && ((NUM_FAILED++))
    ./tools/qa/trapdoor_pylint.py report
    [[ "$?" != 0 ]] && ((NUM_FAILED++))
    ./tools/qa/trapdoor_pep8.py report
    [[ "$?" != 0 ]] && ((NUM_FAILED++))

    git checkout ${CURRENT_BRANCH}

    echo "Number of failed trapdoor checks: ${NUM_FAILED}"
    [[ "${NUM_FAILED}" -gt 0 ]] && exit -1
fi
)
