#!/usr/bin/env bash

# Test if the master branch is present.
echo "Trying to get the commit id of the master branch."
git rev-parse --verify master
RESULT=$?
if [ $RESULT == 0 ]; then
    echo "  OK. Master branch is present."
else
    echo "  NOT OK. Fetching master branch."
    echo "  Assuming it is sufficient to fetch the last 50 commits."
    git config remote.origin.fetch refs/heads/*:refs/remotes/origin/*
    git fetch origin master --depth=50
    git branch master origin/master
fi

# Test if a common ancestor can be found
echo "Trying to get the commit id common ancestor of feature and master."
git merge-base master HEAD
RESULT=$?
if [ $RESULT == 0 ]; then
    echo "  OK. Common ancestor is present."
else
    echo "  NOT OK. Giving up."
    exit -1
fi
