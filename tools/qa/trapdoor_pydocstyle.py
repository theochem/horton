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
"""Trapdoor test using pydocstyle.

This trapdoor uses the pydocstyle program, see http://pydocstyle.readthedocs.io/
The pydocstyle program only tests for a subset of the PEP257, see
https://www.python.org/dev/peps/pep-0257. Not everything can be tested by a program.
"""


from collections import Counter

from trapdoor import TrapdoorProgram, Message, get_source_filenames, run_command


def has_failed(returncode, stdout, stderr):
    """Determine if pydocstyle has failed."""
    return stderr.startswith('Usage:') or stderr.startswith('[Errno')


class PyDocStyleTrapdoorProgram(TrapdoorProgram):
    """A trapdoor program counting the number of pydocstyle messages."""

    def __init__(self):
        """Initialize a PyDocStyleTrapdoorProgram instance."""
        TrapdoorProgram.__init__(self, 'pydocstyle')

    def get_stats(self, config):
        """Run tests using pydocstyle.

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
        # Get version
        command = ['pydocstyle', '--version']
        version = run_command(command, verbose=False)[0].strip()
        print 'USING              : pydocstyle', version

        # Call pydocstyle in the directories containing Python code. All files will be
        # checked, including test files. Missing tests are ignored (D103) because they are
        # already detected by PyLint in a better way.
        command = ['pydocstyle', '--match=.*\\.py', '--add-ignore=D103'] + \
                  config['py_packages'] + \
                  get_source_filenames(config, 'py', unpackaged_only=True)
        output = run_command(command, has_failed=has_failed)[1]

        # Parse the standard output of pydocstyle
        counter = Counter()
        messages = set([])
        lines = output.split('\n')[:-1]
        while len(lines) > 0:
            if 'WARNING: ' in lines[0]:
                lines.pop(0)
            else:
                words = lines.pop(0).split()
                filename, lineno = words[0].split(':')
                code, description = lines.pop(0).split(':', 1)
                code = code.strip()
                description = description.strip()

                key = '%s %s' % (code, filename)
                message = Message(filename, int(lineno), None, '%s %s' % (code, description))

                counter[key] += 1
                messages.add(message)
        return counter, messages


if __name__ == '__main__':
    PyDocStyleTrapdoorProgram().main()
