#!/usr/bin/env python
"""Script that checks if author names were set correctly in commits.

This script ignores and command-line arguments. Just don't provide any. This script
assumes a file AUTHORS is present in the current directory.
"""


import subprocess
import sys


def check_names(key, kind):
    """A generic function to check author/committer names.

    Parameters
    ----------
    key : str
          A key used for formatting the output of git log. Should be 'a' or 'c'.
    kind : str
           Either 'authors' or 'committers'. Should be consistent with the key argument.
    """
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
        print 'UNKNOWN %s:' % kind
        for name in sorted(names):
            print '   ', name
    return len(names) == 0


def main():
    """Check all names against AUTHORS and quit with exit code 1 if an error is found."""
    print 'Checking author and comitter names.'
    pass_a = check_names('a', 'authors')
    pass_c = check_names('c', 'committers')
    if not (pass_a and pass_c):
        sys.exit(1)


if __name__ == '__main__':
    main()
