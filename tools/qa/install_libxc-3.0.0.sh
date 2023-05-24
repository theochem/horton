#!/usr/bin/env bash

VER=3.0.0
NAMEVER=libxc-${VER}
set -e
source tools/qa/common.sh
if [ ! -d "${CACHED}/${NAMEVER}/lib" ]; then
(
    echo -e "${GREEN}Building and installing ${NAMEVER} from scratch${RESET}"
    cd ${QAWORKDIR}
    mkdir -p depbuild
    cd depbuild
    if [ "$VER" = "3.0.0" ]; then
        curl -OL "https://launchpad.net/ubuntu/+archive/primary/+sourcefiles/libxc/3.0.0-1build1/libxc_3.0.0.orig.tar.gz"
    else
        curl -OL "http://www.tddft.org/programs/libxc/down.php?file=${VER}/${NAMEVER}.tar.gz"
    fi
    tar -xzf ${NAMEVER}.tar.gz
    cd ${NAMEVER}
    echo "Actual build and install. This may take a while."
    (CFLAGS='-fPIC' CPPFLAGS='-fPIC' FCCPP=-ffreestanding ./configure --prefix=${CACHED}/${NAMEVER} && make install)
    cp src/funcs_key.c ${CACHED}/${NAMEVER}
#    tail install.log
)
else
    echo -e "${GREEN}Using Cached ${NAMEVER}${RESET}"
fi
