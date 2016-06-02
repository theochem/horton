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
"""Trapdoor test using Cppcheck.

This test calls the cppcheck program, see http://cppcheck.sourceforge.net/.
"""


import subprocess
from xml.etree import ElementTree
from collections import Counter

from trapdoor import TrapdoorProgram, Message, get_source_filenames


class CPPCheckTrapdoorProgram(TrapdoorProgram):
    """A trapdoor program running Cppcheck."""

    def __init__(self):
        """Initialize the CPPCheckTrandoorProgram."""
        TrapdoorProgram.__init__(self, 'cppcheck')

    def get_stats(self, config):
        """Run tests using Cppcheck.

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
        command = ['cppcheck', '--version']
        print 'USING', subprocess.check_output(command, stderr=subprocess.STDOUT).strip()

        # Call Cppcheck
        command = ['cppcheck'] + get_source_filenames(config, 'cpp') + \
                  ['-q', '--enable=all', '--std=c++11', '--xml',
                   '--suppress=missingIncludeSystem', '--suppress=unusedFunction']
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
            text = '%s %s %s' % (error.attrib['severity'], error.attrib['id'], error.attrib['msg'])
            messages.add(Message(error.attrib['file'], int(error.attrib['line']),
                                 None, text))
        return counter, messages


if __name__ == '__main__':
    CPPCheckTrapdoorProgram().main()
