#!/usr/bin/env bash

NAMEVER=$(basename $(dirname "${BASH_SOURCE[0]}"))
set -e
source tools/qa/deps/common.sh
if [ ! -d "${QAWORKDIR}/depinstall/${NAMEVER}/lib/python2.7/site-packages" ]; then
(
    echo -e "${COLOR}Building and installing ${NAMEVER} from scratch${RESET}"
    cd ${QAWORKDIR}
    mkdir -p depbuild
    cd depbuild
    wget http://pypi.python.org/packages/source/h/h5py/h5py-2.2.1.tar.gz
    tar -xzf h5py-2.2.1.tar.gz
    cd h5py-2.2.1
    echo "Actual build and install. This may take a while."
    python setup.py install --prefix=${QAWORKDIR}/depinstall/${NAMEVER} &> install.log
    tail install.log
)
else
    echo -e "${COLOR}Using Cached ${NAMEVER}${RESET}"
fi
