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
"""Trapdoor test using doxygen to find undocumented C++ code.

This test uses doxygen, see http://www.doxygen.org.
"""


import os
import shutil
from collections import Counter

from trapdoor import TrapdoorProgram, Message, run_command


def unwrapped_iter(f):
    """Iterate over unwrapped lines.

    Parameters
    ----------
    f : file
        A text file with wrapped lines. Wrapping is detected by the presence of a
        colon at the end of the line.
    """
    unwrapped_line = ''
    for line in f:
        if line.startswith(' '):
            unwrapped_line += line[:-1]
        else:
            if len(unwrapped_line) > 0:
                yield unwrapped_line
            unwrapped_line = line[:-1]
    if len(unwrapped_line) > 0:
        yield unwrapped_line


class DoxygenTrapdoorProgram(TrapdoorProgram):
    """A trapdoor program counting the number of undocumented C++ functions/methods/..."""

    def __init__(self):
        """Initialize the DoxygenTrapdoorProgram."""
        TrapdoorProgram.__init__(self, 'doxygen')
        self.doxyconf_file = os.path.abspath(os.path.join(self.qaworkdir, 'doxygen.conf'))

    def prepare(self):
        """Make some preparations in feature branch for running doxygen.

        This includes a copy of doc/doxygen.conf to QAWORKDIR.
        """
        TrapdoorProgram.prepare(self)
        shutil.copy('doc/doxygen.conf', self.doxyconf_file)

    def get_stats(self, config):
        """Run tests using doxygen.

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
        command = ['doxygen', '--version']
        print 'USING              : doxygen', run_command(command, verbose=False)[0].strip()

        # Call doxygen in the doc subdirectory, mute output because it only confuses
        command = ['doxygen', self.doxyconf_file]
        run_command(command, cwd=config['doxygen_root'])

        # Parse the file doxygen_warnings log file
        counter = Counter()
        messages = set([])
        prefix = os.getcwd() + '/'

        fn_warnings = os.path.join(config['doxygen_root'], config['doxygen_warnings'])
        with open(fn_warnings, 'r') as f:
            # Doxygen sometimes wraps lines in the warnings log. That is sad, but we
            # have to handle it.
            for line in unwrapped_iter(f):
                location, description = line.split(None, 1)
                filename, lineno = location.split(':')[:2]
                if filename.startswith(prefix):
                    filename = filename[len(prefix):]
                message = Message(filename, int(lineno), None, description)
                counter[filename] += 1
                messages.add(message)
        return counter, messages


if __name__ == '__main__':
    DoxygenTrapdoorProgram().main()
