#!/usr/bin/env bash

if [ "$BUILDKITE_PULL_REQUEST" != "false" ]; then
    source tools/qa/buildkite_common.sh
    checkout_merge_commit
    get_ancestor #writes $SHA variable

    echo "--- Running trapdoors tests"
#    ./tools/qa/check_whitespace.py || report_error "Whitespace errors in some commits"

    TRAPDOORS="trapdoor_coverage.py
    trapdoor_namespace.py"

    for i in ${TRAPDOORS}; do
        .tools/qa/simulate_trapdoor_pr.py -vrA $SHA tools/qa/$i
    done
fi