#!/usr/bin/env bash

set -exo pipefail

abort_error () {
    echo -e "${RED}${1}${RESET}"
    echo -e "${RED}TESTS ABORTED${RESET}"
    exit 1
}

checkout_merge_commit () {
    if [ "$BUILDKITE_PULL_REQUEST" != "false" ]; then
        echo "--- Merging PR"
        API_URL=`echo "$BUILDKITE_PULL_REQUEST_REPO" | sed "s/git:\/\/github.com/https:\/\/api.github.com\/repos/" | sed "s/\.git$/\/pulls\/$BUILDKITE_PULL_REQUEST/"`
        MERGESTATE="null"

        while [ "$MERGESTATE" = "null" ]; do
            MERGESTATE=`curl $API_URL | jq .mergeable`
            if [ "$MERGESTATE" = "true" ]; then
#                git fetch -f origin pull/$BUILDKITE_PULL_REQUEST/merge:temp_merge
                git fetch -f theochem pull/$BUILDKITE_PULL_REQUEST/merge:temp_merge
                git checkout temp_merge
            elif [ "$MERGESTATE" = "null" ]; then
                # git github a chance to compute mergeability
                sleep 5
            else
                abort_error
            fi
        done
    fi

    return 0
}

get_ancestor () {
    if [ "$BUILDKITE_PULL_REQUEST" != "false" ]; then
        echo "--- Finding PR ancestor"
        API_URL=`echo "$BUILDKITE_PULL_REQUEST_REPO" | sed "s/git:\/\/github.com/https:\/\/api.github.com\/repos/" | sed "s/\.git$/\/pulls\/$BUILDKITE_PULL_REQUEST/"`
        SHA=`curl $API_URL | jq .base.sha`
    fi

    return 0
}


# Some colors
GREEN='\e[0;32m'
RED='\e[0;31m'
RESET='\e[0m'

export NUM_FAILED=0

report_error () {
    echo -e "${RED}FAIL: ${1}${RESET}"
    ((NUM_FAILED++))
}
