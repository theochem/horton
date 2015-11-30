#!/usr/bin/env bash

NAMEVER=$(basename $(dirname "${BASH_SOURCE[0]}"))
set -e
source tools/qa/deps/common.sh
if [ ! -d "${QAWORKDIR}/depinstall/${NAMEVER}/lib" ]; then
(
    echo -e "${COLOR}Building and installing ${NAMEVER} from scratch${RESET}"
    cd ${QAWORKDIR}
    mkdir -p depbuild
    cd depbuild
    curl -OL 'http://www.tddft.org/programs/octopus/down.php?file=libxc/libxc-2.2.2.tar.gz'
    tar -xzf libxc-2.2.2.tar.gz
    cd libxc-2.2.2
    echo "Actual build and install. This may take a while."
    (./configure --prefix=${QAWORKDIR}/depinstall/${NAMEVER} --enable-shared && make install) &> install.log
    tail install.log
    echo "Copying libcx_2.2.2/src/funcs_key.c for documentation generator"
    cp src/funcs_key.c ${QAWORKDIR}/depinstall/${NAMEVER}
)
else
    echo -e "${COLOR}Using Cached ${NAMEVER}${RESET}"
fi
