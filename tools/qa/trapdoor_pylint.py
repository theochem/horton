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
'''Trapdoor test using pylint

   This test calls the pep8 program, see http://docs.pylint.org/index.html.
'''


import os
import shutil
import subprocess
from collections import Counter
from trapdoor import TrapdoorProgram


class PylintTrapdoorProgram(TrapdoorProgram):
    def __init__(self):
        TrapdoorProgram.__init__(self, 'pylint')
        self.rcfile = os.path.join(self.qaworkdir, 'pylintrc')

    def initialize(self):
        shutil.copy('tools/qa/.pylintrc', self.rcfile)

    def get_stats(self):
        '''Run tests using Pylint

           Returns
           -------
           counter: collections.Counter
                    counts for different error types in the current checkout
           messages: Set([]) of strings
                     all errors encountered in the current checkout
        '''
        # get Pylint version
        command = ['pylint', '--version', '--rcfile=%s' % self.rcfile]
        version_output = subprocess.check_output(command, stderr=subprocess.STDOUT)
        print 'USING', version_output.split('\n')[0]

        # call Pylint
        command = ['pylint', 'horton', '--rcfile=%s' % self.rcfile]
        command_line = ' '.join(command)
        print 'RUNNING', command_line
        proc = subprocess.Popen(command_line, stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE, shell=True, universal_newlines=True)

        # parse the output of Pylint into standard return values
        lines = proc.stdout.readlines()
        score = lines[-2].split()[6]
        print 'SCORE', score
        counter = Counter()
        messages = set([])
        for line in lines:
            # skip lines that don't contain error messages
            if '.py:' not in line:
                continue
            if line.startswith('Report'):
                break
            # extract error information
            file_name, line = line.strip().split(' ', 1)
            error_id = line[1:].split(']')[0].split(',')[0]
            message = line.split(']')[1].strip()
            counter[error_id] += 1
            messages.add('%-35s  %-40s  %-s' % (error_id, file_name[:-1], message))
        return counter, messages


if __name__ == '__main__':
    PylintTrapdoorProgram().main()
