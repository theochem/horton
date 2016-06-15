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
"""Trapdoor test using nosetests + coverage.

This test calls the nosetests and coverage, see:
* nose: https://nose.readthedocs.org/en/latest/#
* coverage: https://pypi.python.org/pypi/coverage
"""


from collections import Counter
from xml.etree import ElementTree

from trapdoor import TrapdoorProgram, Message, run_command


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
        print 'USING              :', run_command(command, verbose=False)[0].strip()
        command = ['coverage', '--version']
        print 'USING              :', run_command(command, verbose=False)[0].split('\n')[0]

        # Results will be stored in the following variables
        counter = Counter()
        messages = set([])

        # Run fast unit tests with nosetests, with coverage
        command = ['nosetests', '-v', '-a', '!slow', '--with-coverage', '--cover-erase',
                   '--cover-branches',
                   '--cover-package=%s' % ','.join(config['py_packages'])] + \
                   config['py_directories']
        output = run_command(command)[0]
        lines = [line.strip() for line in output.split('\n')]

        # Parse the output of the unit tests
        iline = 0
        for line in lines:
            if len(line) == 0:
                break
            elif line.endswith('FAIL'):
                counter['unit_tests_failed'] += 1
                messages.add(Message(None, None, None, 'nosetests ' + line))
            elif line.endswith('ERROR'):
                counter['unit_tests_error'] += 1
                messages.add(Message(None, None, None, 'nosetests ' + line))
            iline += 1

        # Run the coverage program for a full report. This separate call is needed
        # since coverage-4.1.
        fn_coverage = '%s/coverage.xml' % self.qaworkdir
        command = ['coverage', 'xml', '-o', fn_coverage,
                   '--omit=%s' % ','.join(config['py_test_files'])]
        output = run_command(command)[0]

        # Parse coverage xml output
        et = ElementTree.parse(fn_coverage)
        for class_tag in et.getroot().iter('class'):
            filename = class_tag.attrib['filename']
            for line_tag in class_tag.iter('line'):
                if line_tag.attrib['hits'] == '0':
                    line = int(line_tag.attrib['number'])
                    branch_ends = line_tag.get('missing-branches')
                    if branch_ends is not None:
                        for branch_end in branch_ends.split(','):
                            if branch_end.isdigit():
                                delta = int(branch_end) - line
                                msg = Message(filename, line, None,
                                              'Missed branch to line %+i' % (delta))
                            else:
                                msg = Message(filename, line, None,
                                              'Missed branch to %s' % branch_end)
                            messages.add(msg)
                            counter[filename] += 1
                    messages.add(Message(filename, line, None, 'Missed line'))
                    counter[filename] += 1
        return counter, messages


if __name__ == '__main__':
    CoverageTrapdoorProgram().main()
