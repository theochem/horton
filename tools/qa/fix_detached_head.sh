#!/usr/bin/env bash

# Only do something if the HEAD is detached.
if [ $(git rev-parse --abbrev-ref HEAD) == 'HEAD' ]; then
    # If the current detached HEAD equals the master, then check out the master.
    [[ $(git rev-parse HEAD) == $(git rev-parse master) ]] && git checkout master
fi

# If still not on a branch, give the current branch a name. This must be a feature
# branch or a pull request.
if [ $(git rev-parse --abbrev-ref HEAD) == 'HEAD' ]; then
    # Give the current branch some random name that no serious person would use.
    git checkout -b feature_checkout_2013_foobar
fi
