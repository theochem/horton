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
'''Tool to remove whitespace and tab cruft from source code.'''

import sys


def clean_code(fn):
    '''Return the clean version of a file in which trailing whitespaces on each line are removed,
       tabs are replced by four spcaes, and empty lines at the end of the file are discarded.

       *** Arguments ***

       fn
           The file to be cleaned.
    '''
    print 'Cleaning'.upper(), fn

    # read lines
    with open(fn) as f:
        lines = f.readlines()

    # a flag to be set to true if something really changes. if not, the file
    # is not rewritten.
    changed = False

    # line-by-line stripping of rubish
    for i in xrange(len(lines)):
        # replace tab with four spaces
        line = lines[i].replace('\t', '    ')
        # remove trailing whitespaces from the end of the line
        line = line.rstrip() + '\n'
        # check whether the line has changed
        changed |= (line != lines[i])
        # re-store the line
        lines[i] = line

    # remove empty lines from end of file
    while lines[-1] == '\n':
        changed = True
        lines.pop(-1)

    if changed:
        # re-write lines because the cleaned file is different
        with open(fn, 'w') as f:
            f.writelines(lines)


if __name__ == '__main__':
    # just process all files given as command-line arguments
    for fn in sys.argv[1:]:
        clean_code(fn)
