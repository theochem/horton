#!/usr/bin/env bash

VERSION=1.73
NAMEVER=cppcheck-${VERSION}
set -e
source tools/qa/common.sh
if [ ! -f "${CACHED}/${NAMEVER}/cppcheck" ]; then
(
    echo -e "${GREEN}Building and installing ${NAMEVER} from scratch${RESET}"
    cd ${QAWORKDIR}
    mkdir -p depbuild
    cd depbuild
    curl -OL "http://downloads.sourceforge.net/project/cppcheck/cppcheck/${VERSION}/${NAMEVER}.tar.bz2"
    tar -xjf ${NAMEVER}.tar.bz2
    cd ${NAMEVER}
    echo "Actual build and install. This may take a while."
    make SRCDIR=build CFGDIR=${CACHED}/${NAMEVER}/cfg HAVE_RULES=yes CXXFLAGS="-O2 -DNDEBUG -Wall -Wno-sign-compare -Wno-unused-function" &> install.log
    tail install.log
    mkdir -p ${CACHED}/${NAMEVER}
    cp -av cppcheck ${CACHED}/${NAMEVER}
    cp -av cfg ${CACHED}/${NAMEVER}
)
else
    echo -e "${GREEN}Using Cached ${NAMEVER}${RESET}"
fi
