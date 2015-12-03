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
'''Trapdoor test using cpplint

   This test calls the cpplint program, see
   https://github.com/google/styleguide/tree/gh-pages/cpplint.
'''


import os
import shutil
import stat
import subprocess
from collections import Counter
from glob import glob
from trapdoor import TrapdoorProgram


class CPPLintTrapdoorProgram(TrapdoorProgram):
    def __init__(self):
        TrapdoorProgram.__init__(self, 'cpplint')
        self.cpplint_file = os.path.join(self.qaworkdir, 'cpplint.py')

    def initialize(self):
        shutil.copy('tools/qa/cpplint.py', self.cpplint_file)
        os.chmod(self.cpplint_file, stat.S_IRUSR | stat.S_IWUSR | stat.S_IXUSR)

    def get_stats(self):
        '''Run tests using cpplint

           Returns
           -------
           counter: collections.Counter
                    counts for different error types in the current checkout
           messages: Set([]) of strings
                     all errors encountered in the current checkout
        '''
        # Get version
        print 'Using update #409 of cpplint'

        # Collect all cpp files, except for *_inc.cpp and cext.cpp
        cpp_files = [cpp_file for cpp_file in
            glob('horton/*.cpp') + glob('horton/*.h') +
            glob('horton/*/*.cpp') + glob('horton/*/*.h')
            if not (cpp_file.endswith('_inc.cpp') or cpp_file.endswith('/cext.cpp'))]

        # Call cpplint
        command = [self.cpplint_file] + cpp_files
        print 'RUNNING', ' '.join(command)
        proc = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        output = proc.communicate()[0]

        # Parse the output of cpplint into standard return values
        counter = Counter()
        messages = set([])
        for line in output.split('\n')[:-1]:
            if line.startswith('Done') or line.startswith('Total'):
                continue
            words = line.split()
            filename, lineno = words[0].split(':')[:2]
            description = ' '.join(words[1:-2])
            tag = words[-2]
            priority = words[-1]

            key = '%3s %30s %s' % (priority, tag, filename)
            counter[key] += 1
            messages.add('%3s %30s %40s  %s' % (
                priority,
                tag,
                ('%s:%s' % (filename, lineno)).ljust(40),
                description
            ))

        return counter, messages


if __name__ == '__main__':
    CPPLintTrapdoorProgram().main()
