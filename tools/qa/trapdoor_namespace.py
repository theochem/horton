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
"""Trapdoor test for namespace collisions.

This script imports every module in HORTON and checks for overlapping name spaces
"""


from collections import Counter
import importlib
import os
import sys

from trapdoor import TrapdoorProgram, Message, get_source_filenames


class NamespaceTrapdoorProgram(TrapdoorProgram):
    """A trapdoor program counting namespace collisions."""

    def __init__(self):
        """Initialize the NamespaceTrapdoorProgram."""
        TrapdoorProgram.__init__(self, 'namespace')

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

        # Make sure we test the source tree and not some locally installed copy of HORTON.
        sys.path.insert(0, '.')

        # Find all module names and load namespaces
        namespace = {}
        for filename in get_source_filenames(config, 'py'):
            # Skip irrelevant files
            if filename.endswith('/__init__.py'):
                continue
            if os.path.basename(filename).startswith('test_'):
                continue
            # remove extension and replace / by .
            modulename = filename[:filename.rfind('.')].replace('/', '.')

            # Skip if this modulename is not part of a package
            in_package = False
            for packagename in config['py_packages']:
                if modulename.startswith(packagename):
                    in_package = True
                    break
            if not in_package:
                continue

            # Check all public names of the module.
            module = importlib.import_module(modulename)
            names = dir(module)
            if '__all__' in names:
                for name in names:
                    if name in module.__all__:
                        namespace.setdefault(name, []).append(modulename)
                        if name in config['py_invalid_names']:
                            counter['Invalid name in namespace'] += 1
                            messages.add(Message(filename, None, None,
                                                 'Invalid name in namespace: %s' % name))
            else:
                counter['Missing __all__'] += 1
                messages.add(Message(filename, None, None, 'Missing __all__'))

        # Detect collisions
        for name, modules in namespace.iteritems():
            if len(modules) > 1:
                counter['Namespace collision'] += 1
                text = 'Name \'%s\' found in modules %s' % (name, ' '.join(modules))
                messages.add(Message(None, None, None, text))
        return counter, messages


if __name__ == '__main__':
    NamespaceTrapdoorProgram().main()
