#!/usr/bin/env bash

dir="tools/qa"

if [ "$BUILDKITE_PULL_REQUEST" != "false" ]; then
    source $dir/buildkite_common.sh
    get_ancestor  # Writes $ANCESTOR_SHA variable.

    echo "--- Prep working directory"
    ./cleanfiles.sh

    PATH=$PATH:~/.local/bin  # fix for ubuntu paths

    echo "--- Running trapdoors tests"

    TRAPDOORS="trapdoor_cppcheck.py
    trapdoor_cpplint.py
    trapdoor_doxygen.py
    trapdoor_import.py
    trapdoor_pycodestyle.py
    trapdoor_pydocstyle.py"

    for i in ${TRAPDOORS}; do
        $dir/simulate_trapdoor_pr.py -vA $ANCESTOR_SHA $dir/$i
    done
fi
