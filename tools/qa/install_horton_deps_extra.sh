#!/usr/bin/env bash

# Run all the install scripts in the proper order. After each install, the
# corresponding environment vars are sourced to make sure the following dependencies
# can be built properly.

for DEPDIR in $(cat tools/qa/extra_dependencies.txt); do
(
    tools/qa/${DEPDIR}/install.sh || exit -1
    [[ -f "tools/qa/${DEPDIR}/activate.sh" ]] && source tools/qa/${DEPDIR}/activate.sh
)
done
