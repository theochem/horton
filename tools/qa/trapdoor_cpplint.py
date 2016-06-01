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
"""Trapdoor test using cpplint.

This test calls the cpplint program, see
https://github.com/google/styleguide/tree/gh-pages/cpplint.
"""


import os
import shutil
import stat
import subprocess
from collections import Counter
from fnmatch import fnmatch

from trapdoor import TrapdoorProgram, Message


def get_cpp_files(config):
    """Return a list of .cpp and .h files according to configuration settings.

    This function will search for all .cpp and .h files in the "cpp_directories" and will
    exclude some based in fnmatch patterns provided in the "cpp_exclude" setting.

    Parameters
    ----------
    config : dict
             A dictionary of configuration settings, loaded with json from trapdoor.cfg.
    """
    cpp_files = []
    for cpp_directory in config["cpp_directories"]:
        for dirpath, dirnames, filenames in os.walk(cpp_directory):
            for filename in filenames:
                if (fnmatch(filename, "*.cpp") or fnmatch(filename, "*.h")) and \
                   all(not fnmatch(filename, exclude_filter) for exclude_filter in
                       config["cpp_exclude"]):
                    cpp_files.append(os.path.join(dirpath, filename))
    return cpp_files


class CPPLintTrapdoorProgram(TrapdoorProgram):
    """A trapdoor program running cpplint.py."""

    def __init__(self):
        """Initialize the CPPLintTrapdoorProgram."""
        TrapdoorProgram.__init__(self, 'cpplint')
        self.cpplint_file = os.path.join(self.qaworkdir, 'cpplint.py')

    def prepare(self):
        """Make some preparations in feature branch for running cpplint.py.

        This includes a copy of cpplint.py to QAWORKDIR.
        """
        TrapdoorProgram.prepare(self)
        shutil.copy('tools/qa/cpplint.py', self.cpplint_file)
        os.chmod(self.cpplint_file, stat.S_IRUSR | stat.S_IWUSR | stat.S_IXUSR)

    def get_stats(self, config):
        """Run tests using cpplint.py.

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
        print 'USING cpplint.py update #456'

        # Call cpplint
        command = [self.cpplint_file] + get_cpp_files(config)
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

            counter['%3s %30s %s' % (priority, tag, filename)] += 1
            messages.add(Message(filename, int(lineno), None, '%s %s %s' % (
                priority, tag, description)))

        return counter, messages


if __name__ == '__main__':
    CPPLintTrapdoorProgram().main()
