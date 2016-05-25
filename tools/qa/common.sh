# This file gets sourced by all install and activate scripts

# Some colors
GREEN='\e[0;32m'
RED='\e[0;31m'
RESET='\e[0m'

# Make sure there is a ${QAWORKDIR}
[[ -z ${QAWORKDIR} ]] && export QAWORKDIR=${PWD}/qaworkdir
[[ -d ${QAWORKDIR} ]] || mkdir -p ${QAWORKDIR}

# Fix directories related to caching
export CACHED=${QAWORKDIR}/cached
mkdir -p ${CACHED}

# Make sure matplotlib does not complain about DISPLAY
export MATPLOTLIBRC=${QAWORKDIR}
echo "backend: agg" > $MATPLOTLIBRC/matplotlibrc

# Define an abort function
abort_error() {
    echo -e "${RED}${1}${RESET}"
    echo -e "${RED}TESTS ABORTED (ancestor)${RESET}"
    exit 1
}

# Function for finding and checking ancestor commit
find_ancestor() {
    echo "Getting the common ancestor with the master branch"
    ANCESTOR_COMMIT=$(git merge-base master ${CURRENT_BRANCH})
    local CURRENT_COMMIT=$(git rev-parse HEAD)
    [[ $CURRENT_COMMIT == $ANCESTOR_COMMIT ]] && abort_error "The feature branch is already merged into the master branch."
}
