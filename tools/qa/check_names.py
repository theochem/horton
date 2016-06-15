#!/usr/bin/env python
# -*- coding: utf-8 -*-
# HORTON: Helpful Open-source Research TOol for N-fermion systems.
# Copyright (C) 2011-2016 The HORTON Development Team
#
# This file is part of HORTON.
#
# HORTON is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
#
# HORTON is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, see <http://www.gnu.org/licenses/>
#
# --
"""Script that checks if author names were set correctly in commits.

This script ignores any command-line arguments. Just don't provide any.
This script assumes an AUTHORS file is present in the current directory.
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

    # Turn that list into a set to remove duplicates.
    names = set(names)

    # Make sure the names are in the AUTHORS file.
    with open('AUTHORS', 'r') as f:
        for line in f:
            # chop of the new linefeed character, and
            # remove the author from the names set, if it exist.
            names.discard(line[:-1])
    # Report the names that were not in the AUTHORS file.
    if len(names) != 0:
        print 'UNKNOWN %s:' % kind
        for name in sorted(names):
            print '   ', name
    return len(names) == 0


def main():
    """Check all names against AUTHORS and quit with exit code 1 if an error is found."""
    print 'Checking author and committer names.'.upper()
    pass_a = check_names('a', 'authors')
    pass_c = check_names('c', 'committers')
    if not (pass_a and pass_c):
        sys.exit(1)


if __name__ == '__main__':
    main()
