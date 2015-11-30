# This file gets sourced by all install and activate scripts

# Some colors
COLOR='\033[0;32m'
RESET='\033[0m'

# Make sure there is a ${QAWORKDIR}
[[ -z ${QAWORKDIR} ]] && export QAWORKDIR=${PWD}/qaworkdir
[[ -d ${QAWORKDIR} ]] || mkdir -p ${QAWORKDIR}
