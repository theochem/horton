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
"""Script for removing outdated wheels from cache."""

import os
from glob import glob


# Sort key for the wheels: name and version
def sort_key(s):
    """A sort key for the version in wheel filenames.

    Parameters
    ----------
    s : str
        The wheel filename.
    """
    words = os.path.basename(s).split('-')
    words[1] = words[1].split('.')
    return words


def main():
    """Main program."""
    # Get the directory where we want to remove wheels
    root = os.getenv('PIP_WHEEL_DIR').strip()
    if len(root) == 0:
        raise ValueError('The environment variable PIP_WHEEL_DIR is not defined.')
    if not os.path.isdir(root):
        raise ValueError('The PIP_WHEEL_DIR (%s) does not exist.' % root)
    print 'Removing old wheels in %s' % root

    # get list of all wheels
    all_wheels = sorted(glob('%s/*.whl' % root), key=sort_key)

    # Detect the old wheels, i.e. same name but earlier in the list.
    old_wheels = []
    for i, wheel in enumerate(all_wheels[:-1]):
        name1 = os.path.basename(wheel).split('-')[0]
        name2 = os.path.basename(all_wheels[i+1]).split('-')[0]
        if name1 == name2:
            old_wheels.append(wheel)

    # Remove the old wheels
    for wheel in old_wheels:
        print wheel
        os.remove(wheel)


if __name__ == '__main__':
    main()
