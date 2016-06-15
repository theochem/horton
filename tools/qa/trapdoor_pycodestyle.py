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
"""Trapdoor test using pycodestyle.

This test calls the pycodestyle program, see http://pycodestyle.readthedocs.org/.
This program only covers a part of the PEP8 specification, see
https://www.python.org/dev/peps/pep-0008/. Not everything can be tested by a program.
"""

import os
import shutil
from collections import Counter

import pycodestyle
from trapdoor import TrapdoorProgram, Message


class PyCodeStyleTrapdoorProgram(TrapdoorProgram):
    """A trapdoor program counting the number of pycodestyle messages."""

    def __init__(self):
        """Initialize the PyCodeStyleTrapdoorProgram."""
        TrapdoorProgram.__init__(self, 'pycodestyle')
        self.config_file = os.path.join(self.qaworkdir, 'pycodestyle.ini')

    def prepare(self):
        """Make some preparations in feature branch for running pycodestyle.

        This includes a copy of tools/qa/pycodestyle to QAWORKDIR.
        """
        TrapdoorProgram.prepare(self)
        shutil.copy('tools/qa/%s' % os.path.basename(self.config_file), self.config_file)

    def get_stats(self, config):
        """Run tests using pycodestyle.

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
        print 'USING PYCODESTYLE  :', pycodestyle.__version__

        # Call pycodestyle
        styleguide = pycodestyle.StyleGuide(reporter=CompleteReport,
                                            config_file=self.config_file)
        styleguide.options.exclude.extend(config['py_exclude'])
        print 'EXCLUDED FILES     :', styleguide.options.exclude
        print 'IGNORED MESSAGES   :', styleguide.options.ignore
        print 'MAX LINE LENGTH    :', styleguide.options.max_line_length
        for py_directory in config['py_directories']:
            print 'RUNNING            : pycodestyle %s (through Python API)' % py_directory
            styleguide.input_dir(py_directory)

        # Parse the output of PyCodeStyle into standard return values
        counters = Counter(styleguide.options.report.counters)
        del counters['physical lines']
        del counters['logical lines']
        del counters['files']
        del counters['directories']
        # message on each error
        messages = styleguide.options.report.complete_messages
        assert len(messages) == styleguide.options.report.get_count()
        return counters, messages


class CompleteReport(pycodestyle.StandardReport):
    """Collect and record the results of the checks.

    This subclass is designed to collect all those messages, such that they can be easily
    used in a trapdoor program.
    """

    def __init__(self, options):
        """Initialize a CompleteReport instance."""
        super(CompleteReport, self).__init__(options)
        self.complete_messages = set([])

    def get_file_results(self):
        """Record the result and return the overall count for this file."""
        self._deferred_print.sort()
        for line_number, offset, code, text, _doc in self._deferred_print:
            self.complete_messages.add(Message(
                self.filename, self.line_offset + line_number,
                offset + 1, '%s %s' % (code, text)
            ))


if __name__ == '__main__':
    PyCodeStyleTrapdoorProgram().main()
