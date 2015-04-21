#!/usr/bin/env python
'''Tool to remove whitespace and tab cruft from source code.'''

import sys


def clean_code(fn):
    print 'Cleaning', fn

    # read lines
    with open(fn) as f:
        lines = f.readlines()

    # line-by-line stripping of rubish
    for i in xrange(len(lines)):
        line = lines[i].replace('\t', '    ')
        line = line.rstrip()
        lines[i] = line + '\n'

    # remove empty lines from end of file
    while lines[-1] == '\n':
        lines.pop(-1)

    # write lines
    with open(fn, 'w') as f:
        f.writelines(lines)


if __name__ == '__main__':
    # just process all files given as command-line arguments
    for fn in sys.argv[1:]:
        clean_code(fn)
