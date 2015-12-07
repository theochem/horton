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
'''Trapdoor test for PEP257

See https://www.python.org/dev/peps/pep-0257 for the complete PEP.

This script script uses the pep257 package which only tests for a subset of the PEP257
guidelines.
'''


import os
import subprocess
from collections import Counter
from trapdoor import TrapdoorProgram


class PEP257TrapdoorProgram(TrapdoorProgram):
    def __init__(self):
        TrapdoorProgram.__init__(self, 'pep257')

    def get_stats(self, config):
        '''Run tests using pep257

           Returns
           -------
           counter: collections.Counter
                    counts for different error types in the current checkout
           messages: Set([]) of strings
                     all errors encountered in the current checkout
        '''
        # Get version
        command = ['pep257', '--version']
        print 'USING pep257', subprocess.check_output(command).strip()

        # Call doxygen in the doc subdirectory
        command = ['pep257'] + config['py_directories']
        print 'RUNNING', ' '.join(command)
        proc = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        output = proc.communicate()[0]

        # Parse the file doc/doxygen_warnings.log
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
                message = '%s  %-40s  %s' % (
                    code,
                    '%s:%s' % (filename, lineno),
                    description
                )

                counter[key] += 1
                messages.add(message)
        return counter, messages


if __name__ == '__main__':
    PEP257TrapdoorProgram().main()
