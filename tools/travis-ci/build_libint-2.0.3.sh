#!/usr/bin/env bash

# exit script when a line fails (instead of going to the next line)
set -e

COLOR='\033[0;32m'
RESET='\033[0m'

# check to see if libint-2.0.3 is already present
if [ ! -d "$HOME/depinstall/libint-2.0.3/lib" ]; then
(
    echo -e "${COLOR}Building and installing LIBINT 2.0.3 from scratch${RESET}"
    cd
    mkdir -p depbuild
    cd depbuild
    wget 'http://downloads.sourceforge.net/project/libint/libint-for-mpqc/libint-2.0.3-stable.tgz'
    tar -xzf libint-2.0.3-stable.tgz
    cd libint-2.0.3-stable
    echo "Actual build and install. This may take a while."
    echo "Every 100th line of build output is shown to keep Travis-CI happy without generating too much output"
    (CONFIG_SHELL=$(which sh) ./configure --prefix=${HOME}/depinstall/libint-2.0.3 --enable-shared && make install) | awk '!(NR % 100)'
)
else
    echo -e "${COLOR}Using Cached LIBINT 2.0.3${RESET}"
fi

find $HOME/depinstall/libint-2.0.3
