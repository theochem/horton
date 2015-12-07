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
#--
'''Trapdoor test using pep8

   This test calls the pep8 program, see http://pep8.readthedocs.org/.
'''

import os
import pep8
import shutil
from collections import Counter
from trapdoor import TrapdoorProgram


class PEP8TrapdoorProgram(TrapdoorProgram):
    def __init__(self):
        TrapdoorProgram.__init__(self, 'pep8')
        self.config_file = os.path.join(self.qaworkdir, 'pep8')

    def initialize(self):
        TrapdoorProgram.initialize(self)
        shutil.copy('tools/qa/pep8', self.config_file)

    def get_stats(self, config):
        '''Run tests using pep8

           Returns
           -------
           counter: collections.Counter
                    counts for different error types in the current checkout
           messages: Set([]) of strings
                     all errors encountered in the current checkout
        '''
        # Get version
        print 'USING pep8', pep8.__version__

        # Call Pep8
        pep8check = pep8.StyleGuide(reporter=CompleteReport, config_file=self.config_file)
        print 'EXCLUDED FILES', pep8check.options.exclude
        print 'IGNORED MESSAGES', pep8check.options.ignore
        print 'MAX LINE LENGTH', pep8check.options.max_line_length
        for py_directory in config['py_directories']:
            print 'RUNNING pep8', py_directory
            pep8check.input_dir(py_directory)

        # Parse the output of Pep8 into standard return values
        counters = Counter(pep8check.options.report.counters)
        del counters['physical lines']
        del counters['logical lines']
        del counters['files']
        del counters['directories']
        # message on each error
        messages = set(pep8check.options.report.complete_message)
        assert len(messages) == pep8check.options.report.get_count()
        return counters, messages


class CompleteReport(pep8.StandardReport):
    '''
    Collect and record the results of the checks.
    '''
    def __init__(self, options):
        super(CompleteReport, self).__init__(options)
        self.complete_message = []

    def get_file_results(self):
        '''Record the result and return the overall count for this file.'''
        self._deferred_print.sort()
        for line_number, offset, code, text, doc in self._deferred_print:
            # record the error message specifications.
            message = '%15s  %50s  %s' % (
                code,
                ('%s:%s:%s' % (self.filename, self.line_offset + line_number, offset + 1)).ljust(50),
                text)
            if self._show_source:
                if line_number > len(self.lines):
                    line = ''
                else:
                    line = self.lines[line_number - 1]
            self.complete_message.append(message)


if __name__ == '__main__':
    PEP8TrapdoorProgram().main()
