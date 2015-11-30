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
    wget 'http://downloads.sourceforge.net/project/libint/libint-for-mpqc/libint-2.0.3-stable.tgz'
    tar -xzf libint-2.0.3-stable.tgz
    cd libint-2.0.3-stable
    echo "Actual build and install. This may take a while."
    echo "Every 100th line of the build output is shown to keep Travis-CI happy without generating too much output"
    (CONFIG_SHELL=$(which sh) ./configure --prefix=${QAWORKDIR}/depinstall/${NAMEVER} --enable-shared && make install) 2>&1 | tee install.log | awk '!(NR % 100)'
    tail install.log
)
else
    echo -e "${COLOR}Using Cached ${NAMEVER}${RESET}"
fi
