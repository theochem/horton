#!/usr/bin/env bash

NAMEVER=$(basename $(dirname "${BASH_SOURCE[0]}"))
set -e
source tools/qa/common.sh
if [ ! -d "${CACHED}/${NAMEVER}/lib" ]; then
(
    echo -e "${GREEN}Building and installing ${NAMEVER} from scratch${RESET}"
    cd ${QAWORKDIR}
    mkdir -p depbuild
    cd depbuild
    curl -OL 'http://www.tddft.org/programs/octopus/down.php?file=libxc/libxc-2.2.2.tar.gz'
    tar -xzf libxc-2.2.2.tar.gz
    cd libxc-2.2.2
    echo "Actual build and install. This may take a while."
    (./configure --prefix=${CACHED}/${NAMEVER} --enable-shared && make install) &> install.log
    tail install.log
    echo "Copying libcx_2.2.2/src/funcs_key.c for documentation generator"
    cp src/funcs_key.c ${CACHED}/${NAMEVER}
)
else
    echo -e "${GREEN}Using Cached ${NAMEVER}${RESET}"
fi
