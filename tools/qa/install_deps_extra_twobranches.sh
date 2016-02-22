#!/usr/bin/env bash

# 1) Build all dependencies of the current branch. This is always needed. If the build
#    fails, the script exits with -1, letting the build bot know that it should stop.

tools/qa/install_deps_extra.py || exit -1

# 2) Get the current branch name. If not master, then check out master and build all
#    dependencies. Finally check out the current branch again.

CURRENT_BRANCH=$(git rev-parse --abbrev-ref HEAD)
if [ "${CURRENT_BRANCH}" != 'master' ]; then
    git checkout master
    tools/qa/install_deps_extra.py || exit -1
    git checkout ${CURRENT_BRANCH}
fi

# Each install script called above must follow some conventions:
#
# - It installs the dependency in ${QAWORKDIR}/depinstall/name-version/
#
# - It only builds the dependency if not present yet
#
# - A script activate_name_version.sh is provided that can be sourced before the
#   dependency can be used. This activate script is called before running the actual
#   tests. It should not do anything if the corresponding depinstall directory is not
#   present, which is useful when a developer has manually installed the dependency
#   elsewhere.
