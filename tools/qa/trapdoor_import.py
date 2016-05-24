#!/usr/bin/env python
# -*- coding: utf-8 -*-
# HORTON: Helpful Open-source Research TOol for N-fermion systems.
# Copyright (C) 2011-2015 The HORTON Development Team
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
"""Trapdoor test for import conventions.

This script counts the number of bad imports. The following is not allowed in a package:

* Importing from the top-level of its own package, e.g.:

  .. code-block:: python

        from package import foo
"""


import os
import codecs
from collections import Counter

from trapdoor import TrapdoorProgram


class ImportTrapdoorProgram(TrapdoorProgram):
    """A trapdoor program counting the number of bad imports."""

    def __init__(self):
        """Initialize the ImportTrapdoorProgram."""
        TrapdoorProgram.__init__(self, 'import')

    def get_stats(self, config):
        """Count number of bad tests.

        Parameters
        ----------
        config : dict
                 The dictionary loaded from ``trapdoor.cfg``.

        Returns
        -------
        counter : collections.Counter
                  Counts of the number of messages of a specific type in a certain file.
        messages : Set([]) of strings
                   All errors encountered in the current branch.
        """
        # Output variables
        counter = Counter()
        messages = set([])

        # Loop all python and cython files
        for pydir in config['py_directories']:
            for rootdir, subdirs, filenames in os.walk(pydir):
                for filename in filenames:
                    # Only consider relevant files
                    if not (filename.endswith('.py') or filename.endswith('.pyx')):
                        continue
                    if filename.startswith('test_') or filename == '__init__.py':
                        continue
                    # Look for bad imports
                    path = os.path.join(rootdir, filename)
                    with codecs.open(path, encoding='utf-8') as f:
                        for lineno, line in enumerate(f):
                            for py_package in config['py_packages']:
                                if u'from %s import' % py_package in line:
                                    counter['Wrong imports in %s' % path] += 1
                                    messages.add('Wrong import: %s:%i' % (path, lineno))

        return counter, messages


if __name__ == '__main__':
    ImportTrapdoorProgram().main()
