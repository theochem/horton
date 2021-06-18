#!/usr/bin/env bash

http://www.tddft.org/programs/libxc/down.php?file=2.2.2/libxc-2.2.2.tar.gz

VER=2.2.2
NAMEVER=libxc-${VER}
set -e
source tools/qa/common.sh
if [ ! -d "${CACHED}/${NAMEVER}/lib" ]; then
(
    echo -e "${GREEN}Building and installing ${NAMEVER} from scratch${RESET}"
    cd ${QAWORKDIR}
    mkdir -p depbuild
    cd depbuild
    curl -OL "http://www.tddft.org/programs/libxc/down.php?file=${VER}/${NAMEVER}.tar.gz"
    tar -xzf libxc-2.2.2.tar.gz
    cd libxc-2.2.2
    echo "Actual build and install. This may take a while."
    (CFLAGS='-fPIC' CPPFLAGS='-fPIC' FCCPP=-ffreestanding ./configure --prefix=${CACHED}/${NAMEVER} && make install) &> install.log
    cp src/funcs_key.c ${CACHED}/${NAMEVER}
    tail install.log
)
else
    echo -e "${GREEN}Using Cached ${NAMEVER}${RESET}"
fi
