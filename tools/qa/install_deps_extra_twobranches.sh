#!/usr/bin/env bash

source tools/qa/common.sh

# 1) Build all dependencies of the current branch. This is always needed. If the build
#    fails, the script exits with -1, letting the build bot know that it should stop.

tools/qa/install_deps_extra.py || exit -1

# 2) Get the current branch name. If not master, then check out the common ancestor with
#    master and build all dependencies. Finally check out the current branch again.

CURRENT_BRANCH=$(git rev-parse --abbrev-ref HEAD)
ANCESTOR_COMMIT=$(git merge-base master ${CURRENT_BRANCH})
CURRENT_COMMIT=$(git rev-parse HEAD)
if [ "${CURRENT_BRANCH}" != 'master' ] && [ "${CURRENT_COMMIT}" != ${ANCESTOR_COMMIT} ]; then
        git checkout ${ANCESTOR_COMMIT}
        tools/qa/install_deps_extra.py || exit -1
        git checkout ${CURRENT_BRANCH}
else
    echo "No need to install dependencies of ancestor because HEAD is (merged in) master."
fi

# Each install script called above must follow some conventions:
#
# - It installs the dependency in ${QAWORKDIR}/depinstall/name-version/
#
# - It only builds the dependency if not present yet
