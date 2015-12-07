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
'''Trapdoor test using Cppcheck

   This test calls the cppcheck program, see http://cppcheck.sourceforge.net/.
'''


import subprocess
from xml.etree import ElementTree
from collections import Counter
from trapdoor import TrapdoorProgram
from trapdoor_cpplint import get_cpp_files


class CPPCheckTrapdoorProgram(TrapdoorProgram):
    def __init__(self):
        TrapdoorProgram.__init__(self, 'cppcheck')

    def get_stats(self, config):
        '''Run tests using Cppcheck

           Returns
           -------
           counter: collections.Counter
                    counts for different error types in the current checkout
           messages: Set([]) of strings
                     all errors encountered in the current checkout
        '''
        # Get version
        command = ['cppcheck', '--version']
        print 'USING', subprocess.check_output(command, stderr=subprocess.STDOUT).strip()

        # Call Cppcheck
        cpp_files = get_cpp_files(config)
        command = ['cppcheck'] + cpp_files + ['-q', '--enable=all',
                   '--std=c++11', '--xml', '--suppress=missingIncludeSystem']
        print 'RUNNING', ' '.join(command)
        xml_str = subprocess.check_output(command, stderr=subprocess.STDOUT)
        etree = ElementTree.fromstring(xml_str)

        # Parse the output of Cppcheck into standard return values
        counter = Counter()
        messages = set([])
        for error in etree:
            if 'file' not in error.attrib:
                continue
            key = '%15s  %40s  %30s' % (
                error.attrib['severity'],
                error.attrib['file'].ljust(40),
                error.attrib['id'].ljust(30),
            )
            counter[key] += 1
            messages.add('%15s  %40s  %s' % (
                error.attrib['severity'],
                ('%s:%s' % (error.attrib['file'], error.attrib['line'])).ljust(40),
                error.attrib['msg'],
            ))
        return counter, messages


if __name__ == '__main__':
    CPPCheckTrapdoorProgram().main()
