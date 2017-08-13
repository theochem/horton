#!/usr/bin/env bash

NAMEVER=libxc-3.0.0
set -e
source tools/qa/common.sh
if [ ! -d "${CACHED}/${NAMEVER}/lib" ]; then
(
    echo -e "${GREEN}Building and installing ${NAMEVER} from scratch${RESET}"
    cd ${QAWORKDIR}
    mkdir -p depbuild
    cd depbuild
    curl -OL "http://www.tddft.org/programs/octopus/down.php?file=libxc/${NAMEVER}.tar.gz"
    tar -xzf ${NAMEVER}.tar.gz
    cd ${NAMEVER}
    echo "Actual build and install. This may take a while."
    (CFLAGS='-fPIC' CPPFLAGS='-fPIC' FCCPP=-ffreestanding ./configure --prefix=${CACHED}/${NAMEVER} && make install) &> install.log
    cp src/funcs_key.c ${CACHED}/${NAMEVER}
    tail install.log
)
else
    echo -e "${GREEN}Using Cached ${NAMEVER}${RESET}"
fi
