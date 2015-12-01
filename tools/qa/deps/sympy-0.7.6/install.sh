#!/usr/bin/env bash

NAMEVER=$(basename $(dirname "${BASH_SOURCE[0]}"))
set -e
source tools/qa/deps/common.sh
if [ ! -d "${QAWORKDIR}/depinstall/${NAMEVER}/lib/python2.7/site-packages" ]; then
(
    echo -e "${GREEN}Building and installing ${NAMEVER} from scratch${RESET}"
    cd ${QAWORKDIR}
    mkdir -p depbuild
    cd depbuild
    wget http://pypi.python.org/packages/source/s/sympy/sympy-0.7.6.tar.gz
    tar -xzf sympy-0.7.6.tar.gz
    cd sympy-0.7.6
    echo "Actual build and install. This may take a while."
    python setup.py install --prefix=${QAWORKDIR}/depinstall/${NAMEVER} &> install.log
    tail install.log
)
else
    echo -e "${GREEN}Using Cached ${NAMEVER}${RESET}"
fi
