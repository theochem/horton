#!/usr/bin/env python
'''Script that checks if author names were set correctly in commits'''


import subprocess
import sys


def check_names(key, kind):
    # Get a list of authors/committer names for every commit.
    command = ['git', 'log', '--format=%%%sN <%%%sE>' % (key, key)]
    names = subprocess.check_output(command).split('\n')[:-1]

    # Turn that list into a set
    names = set(names)

    # Make sure the names are in the AUTHORS file
    with open('AUTHORS', 'r') as f:
        for line in f:
            # chop of the new-line
            names.discard(line[:-1])
    if len(names) != 0:
        print 'Unknown %s:' % kind
        for name in sorted(names):
            print '   ', name
    return len(names) == 0


def main():
    pass_a = check_names('a', 'authors')
    pass_c = check_names('c', 'committers')
    if not (pass_a and pass_c):
        sys.exit(1)


if __name__ == '__main__':
    main()
