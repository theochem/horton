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
"""Trapdoor test using nosetests + coverage.

This test calls the nosetests and coverage, see:
* nose: https://nose.readthedocs.org/en/latest/#
* coverage: https://pypi.python.org/pypi/coverage
"""


import subprocess
from collections import Counter
from trapdoor import TrapdoorProgram


class CoverageTrapdoorProgram(TrapdoorProgram):
    """A trapdoor program running nosetests with coverage analysis."""

    def __init__(self):
        """Initialize the CoverageTrapdoorProgram."""
        TrapdoorProgram.__init__(self, 'coverage')

    def get_stats(self, config):
        """Run tests using nosetests with coverage analysis.

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
        command = ['nosetests', '--version']
        print 'USING', subprocess.check_output(command, stderr=subprocess.STDOUT).strip()
        command = ['coverage', '--version']
        print 'USING', subprocess.check_output(command, stderr=subprocess.STDOUT).split('\n')[0]

        # Run fast unit tests with nosetests, with coverage
        command = ['nosetests', '-v', '-a', '!slow', '--with-coverage', '--cover-erase',
                   '--cover-branches',
                   '--cover-package=%s' % ','.join(config['py_packages'])] + \
                   config['py_directories']
        print 'RUNNING', ' '.join(command)
        proc = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        output = proc.communicate()[0]
        lines = [line.strip() for line in output.split('\n')]

        # Results will be stored in the following variables
        counter = Counter()
        messages = set([])

        # Parse the output of the unit tests
        iline = 0
        for line in lines:
            if len(line) == 0:
                break
            elif line.endswith('FAIL'):
                counter['unit_tests_failed'] += 1
                messages.add('nosetests ' + line)
            elif line.endswith('ERROR'):
                counter['unit_tests_error'] += 1
                messages.add('nosetests ' + line)
            iline += 1

        # Run the coverage program for a full report. This separate call is needed
        # since coverage-4.1.
        command = ['coverage', 'report', '-m']
        print 'RUNNING', ' '.join(command)
        proc = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        output = proc.communicate()[0]
        lines = [line.strip() for line in output.split('\n')]

        # Parse coverage output
        # Remove the first two lines
        lines = lines[2:]
        # Process line by line
        for line in lines:
            if line.startswith('--------'):
                break
            words = line.split()
            miss = int(words[2])  # number of lines missed
            if miss > 0:
                filename = words[0]
                counter['missed lines in ' + filename] += miss
                for linenos in words[4:]:
                    linenos = linenos.strip()
                    if linenos.endswith(','):
                        linenos = linenos[:-1]
                    if len(linenos) > 0:
                        messages.add('coverage %s %s' % (filename, linenos))

        return counter, messages


if __name__ == '__main__':
    CoverageTrapdoorProgram().main()
