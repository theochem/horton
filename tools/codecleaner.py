#!/usr/bin/env python
'''Tool to remove whitespace and tab cruft from source code.'''

import sys


def clean_code(fn):
    print 'Cleaning', fn

    # read lines
    with open(fn) as f:
        lines = f.readlines()

    # this will be set to true if something really changes. if not, the file
    # is not rewritten.
    changed = False

    # line-by-line stripping of rubish
    for i in xrange(len(lines)):
        line = lines[i].replace('\t', '    ')
        line = line.rstrip() + '\n'
        changed |= (line != lines[i])
        lines[i] = line

    # remove empty lines from end of file
    while lines[-1] == '\n':
        changed = True
        lines.pop(-1)

    if changed:
        # write lines
        with open(fn, 'w') as f:
            f.writelines(lines)


if __name__ == '__main__':
    # just process all files given as command-line arguments
    for fn in sys.argv[1:]:
        clean_code(fn)
