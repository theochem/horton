#!/usr/bin/env bash

# exit script when a line fails (instead of going to the next line)
set -e

COLOR='\033[0;32m'
RESET='\033[0m'

# check to see if libxc-2.2.2 is already present
if [ ! -d "$HOME/depinstall/libxc-2.2.2/lib" ]; then
(
    echo -e "${COLOR}Building and installing LibXC 2.2.2 from scratch${RESET}"
    cd
    mkdir -p depbuild
    cd depbuild
    curl -OL 'http://www.tddft.org/programs/octopus/down.php?file=libxc/libxc-2.2.2.tar.gz'
    tar -xzf libxc-2.2.2.tar.gz
    cd libxc-2.2.2
    echo "Actual build and install. This may take a while."
    (./configure --prefix=${HOME}/depinstall/libxc-2.2.2 --enable-shared && make install) &> install.log
    tail install.log
)
else
    echo -e "${COLOR}Using Cached LibXC 2.2.2${RESET}"
fi

find $HOME/depinstall/libxc-2.2.2
