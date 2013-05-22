#!/usr/bin/env python

# This scripts separates all words with a fixed spacing of 15 characters, if the
# line contains more than two words. This is used to format the file data/elements.txt

import sys

with open(sys.argv[1]) as f:
    lines = f.readlines()

for i in xrange(len(lines)):
    words = lines[i].split()
    if len(words) > 2 and not lines[i].strip().startswith('#'):
        lines[i] = ' '.join('%14s' % word for word in words) + '\n'

with open(sys.argv[1], 'w') as f:
    f.writelines(lines)
