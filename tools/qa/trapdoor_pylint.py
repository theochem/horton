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
"""Trapdoor test using pylint.

This test calls the pylint program, see http://docs.pylint.org/index.html.
"""


import os
import shutil
import subprocess
from collections import Counter
from trapdoor import TrapdoorProgram


class PylintTrapdoorProgram(TrapdoorProgram):
    """A trapdoor program counting the number of pylint messages."""

    def __init__(self):
        """Initialize a PylintTrapdoorProgram instance."""
        TrapdoorProgram.__init__(self, 'pylint')
        self.rcfile = os.path.join(self.qaworkdir, 'pylintrc')

    def prepare(self):
        """Make some preparations in feature branch for running pylint.

        This includes a copy of tools/qa/pylintrc to QAWORKDIR.
        """
        TrapdoorProgram.prepare(self)
        shutil.copy('tools/qa/pylintrc', self.rcfile)

    def get_stats(self, config):
        """Run tests using Pylint.

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
        # get Pylint version
        command = ['pylint', '--version', '--rcfile=%s' % self.rcfile]
        version_output = subprocess.check_output(command, stderr=subprocess.STDOUT)
        print 'USING   :', ''.join(version_output.split('\n')[:2])

        # call Pylint
        command = ['pylint'] + config['py_directories'] + ['--rcfile=%s' % self.rcfile]
        print 'RUNNING :', ' '.join(command)
        proc = subprocess.Popen(command, stdout=subprocess.PIPE)

        # parse the output of Pylint into standard return values
        output = proc.communicate()[0]
        lines = output.split('\n')[:-1]
        score = lines[-2].split()[6]
        print 'SCORE   :', score
        counter = Counter()
        messages = set([])
        for line in lines:
            # skip lines that don't contain error messages
            if '.py:' not in line:
                continue
            if line.startswith('Report'):
                break
            # extract error information
            msg_id, location, msg = line.split(' ', 2)
            counter[msg_id] += 1
            messages.add('%-5s  %-40s  %s' % (msg_id, location, msg))
        return counter, messages


if __name__ == '__main__':
    PylintTrapdoorProgram().main()
