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


abort_error() {
    echo -e "${RED}${1}${RESET}"
    echo -e "${RED}TESTS ABORTED${RESET}"
    exit 1
}
