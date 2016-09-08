#!/usr/bin/env bash

if [ "$BUILDKITE_PULL_REQUEST" != "false" ]; then
    source tools/qa/buildkite_common.sh
    checkout_merge_commit
    get_ancestor #writes $SHA variable

    echo "--- Build refatoms"
    rm -rf data/refatoms/*.h5 data/refatoms/*.tar.bz2
    make -C data/refatoms/

    echo "--- Running trapdoors tests"
    ./tools/qa/check_whitespace.py || report_error "Whitespace errors in some commits"

    TRAPDOORS="trapdoor_coverage.py
    trapdoor_namespace.py"

    export PYTHONPATH=`pwd`
    export HORTONDATA=`pwd`/data

    for i in ${TRAPDOORS}; do
        tools/qa/simulate_trapdoor_pr.py -vrA $SHA tools/qa/$i
    done
fi
