#!/usr/bin/env bash

NAMEVER=libint-2.0.3
set -e
source tools/qa/common.sh
if [ ! -d "${CACHED}/${NAMEVER}/lib" ]; then
(
    echo -e "${GREEN}Building and installing ${NAMEVER} from scratch${RESET}"
    cd ${QAWORKDIR}
    mkdir -p depbuild
    cd depbuild
    curl -OL 'http://downloads.sourceforge.net/project/libint/libint-for-mpqc/libint-2.0.3-stable.tgz'
    tar -xzf libint-2.0.3-stable.tgz
    cd libint-2.0.3-stable
    echo "Actual build and install. This may take a while."
    echo "Every 100th line of the build output is shown to keep Travis-CI happy without generating too much output"
    (CONFIG_SHELL=$(which sh) CFLAGS='-fPIC' CPPFLAGS='-fPIC' ./configure --with-cxx-optflags='-O1' --prefix=${CACHED}/${NAMEVER} && make install) 2>&1 | tee install.log | awk '!(NR % 100)'
    tail install.log
)
else
    echo -e "${GREEN}Using Cached ${NAMEVER}${RESET}"
fi
