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
'''Trapdoor test using nosetests + coverage

   This test calls the nosetests and coverage, see
       nose: https://nose.readthedocs.org/en/latest/#
       coverage: https://pypi.python.org/pypi/coverage
'''


import subprocess
from collections import Counter
from trapdoor import main


def get_stats_coverage_check():
    '''Run tests using Cppcheck

       Returns
       -------
       counter: collections.Counter, {string: float}
                string: the name of nosetests
                float: percentage of coverage
       messages: Set([]) of strings
                name of nosetests: [uncovered lines]
                 
    '''
    # Get version
    command = ['nosetests', '--version']
    print 'Using', subprocess.check_output(command, stderr=subprocess.STDOUT).strip()
    trial_command = ['coverage', 'coverage-2.7', 'coverage27']
    command = []
    for i in trial_command:
        try:
            subprocess.check_output(i)
        except:
            continue
        else:
            cover_command = i
            command = [cover_command, '--version']
    if command:
        print 'Using', subprocess.check_output(command, stderr=subprocess.STDOUT).strip()
    else:
        print("excutable coverage command not found")
        return
    # Call coverage
    command = ['nosetests', '--with-coverage', '--cover-erase',
               '--cover-package=horton.io', 'horton.io']
    print 'RUNNING', ' '.join(command)
    result_str = subprocess.check_output(command, stderr=subprocess.STDOUT)
    print result_str

    # Parse the output of Cppcheck into standard return values
    counter = Counter()
    messages = set([])
    lines_result = result_str.strip(".").strip().split("\n")
    for eachline in lines_result[2:-4]:
        line_elements = eachline.replace(", ", ",").split()
        # print line_elements
        if len(line_elements) == 5:
            key = line_elements[0]
            value = float(line_elements[3].strip("%")) / 100
            info = line_elements[4]
            counter[key] = value
            messages.add("{0}: {1}".format(key, info))
    return counter, messages


if __name__ == '__main__':
    main(get_stats_coverage_check)
