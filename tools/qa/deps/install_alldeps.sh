#!/usr/bin/env bash

# Run all the install scripts in the proper order. After each install, the
# corresponding environment vars are sourced to make sure the following dependencies
# can be built properly.

for DEPDIR in $(cat tools/qa/deps/dirs.txt); do
(
    tools/qa/deps/${DEPDIR}/install.sh || exit -1
    [[ -f "tools/qa/deps/${DEPDIR}/activate.sh" ]] && source tools/qa/deps/${DEPDIR}/activate.sh
)
done
