#!/usr/bin/env bash

# Run two sets of tests and exit nonzero if one fails.
./tools/qa/test_all_feature.sh || exit 1
./tools/qa/test_all_ancestor.sh || exit 1
