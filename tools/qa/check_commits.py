#!/usr/bin/env python
# -*- coding: utf-8 -*-
# HORTON: Helpful Open-source Research TOol for N-fermion systems.
# Copyright (C) 2011-2017 The HORTON Development Team
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
"""Check commits for bad habits. Can be used as pre-commit when no arguments are given."""

import argparse
import subprocess
import sys


def main():
    """Main program."""
    ancestor = parse_args()
    check_commits(ancestor)


def parse_args():
    """Determine which commits to check, from command-line args.

    When no arguments are given, the index (commit to be made) will be checked. Otherwise,
    all commits up to the ancestor are checked.
    """
    parser = argparse.ArgumentParser(description='Check basic formatting conventions in commits.')
    parser.add_argument('ancestor', type=str, default=None, nargs='?',
                        help='The ancestor up to which to check the commits. This can be '
                             'an SHA1 hash of the commit or the name of a branch. In '
                             'the latter case the HEAD of that branch is used as '
                             'ancestor. When not given, the index (containing files to '
                             'be committed) is checked.')
    args = parser.parse_args()
    return args.ancestor


def check_commits(ancestor):
    """Check commits up to ancestor for bad habits.

    When ancestor is None, the staged files are checked.
    """
    # If all goes well, the exit code remains zero.
    result = 0

    if ancestor is None:
        print 'Checking for bad habits in files to be committed ...'
        diff_args = [None]
    else:
        print 'Checking for bad habits in commits up to {} ...'.format(ancestor)
        lines = subprocess.check_output(
            ['git', 'log', '%s..HEAD' % ancestor, '--pretty=oneline', '--color=never'],
            stderr=subprocess.STDOUT).decode("utf-8").splitlines()

        # Parse the output
        commit_ids = []
        for line in lines:
            commit_id = line.partition(' ')[0]
            commit_ids.append(commit_id)
        commit_ids.append(ancestor)

        # Define the ranges to be checked with git diff
        diff_args = []
        for i in xrange(len(commit_ids) - 1):
            diff_args.append('{}..{}'.format(commit_ids[i + 1], commit_ids[i]))

    # Loop over all commits and check them:
    for diff_arg in diff_args:
        if ancestor is not None:
            print '{}'.format(diff_arg)
        # Get a list of files in the diff.
        status_lines = subprocess.check_output(
            ['git', 'diff', '--name-status', '-M', diff_arg or '--cached'],
            stderr=subprocess.STDOUT).decode("utf-8").splitlines()

        status_lines = [line for line in status_lines if
                        not (line.endswith(".npy") or line.endswith(".npz"))]

        # Check all files in the diff.
        for status_line in status_lines:
            # Parse the status line
            words = status_line.split()
            status = words[0][0]
            if status == 'D':
                continue
            elif status == 'R':
                old_filename, new_filename = words[1:]
            else:
                new_filename = words[1]
                old_filename = new_filename

            # Get the old and new blobs
            if diff_arg is None:
                old_commit_id = 'HEAD'
                new_commit_id = ''
            else:
                old_commit_id, new_commit_id = diff_arg.split('..')
            if status == 'A':
                old_blob = '--'
            else:
                old_blob = '{}:{}'.format(old_commit_id, old_filename)
            new_blob = '{}:{}'.format(new_commit_id, new_filename)

            # Get the diff
            if status == 'R':
                print '   {} {} -> {}'.format(status, old_filename, new_filename)
            else:
                print '   {} {}'.format(status, new_filename)
            diff_lines = subprocess.check_output(
                ['git', 'diff', '--color=never', '-M', old_blob, new_blob],
                stderr=subprocess.STDOUT).decode("utf-8").splitlines()

            # Run some tests on the diff.
            for line in diff_lines:
                if line.startswith('diff') or line.startswith('index') or \
                        line.startswith('-') or line.startswith('+++'):
                    pass
                elif line.startswith('@@'):
                    line_number = line.split()[2]
                    if ',' in line_number:
                        line_number = line_number.partition(',')[0]
                    line_number = int(line_number)
                elif line.startswith('+'):
                    location = '{}:{}'.format(new_filename, line_number)
                    if '\t' in line and not new_filename.endswith('Makefile'):
                        print '      Tab                   {}'.format(location)
                        result = 1
                    if line.endswith(' '):
                        print '      Trailing whitespace   {}'.format(location)
                        result = 1
                    line_number += 1
                elif line.startswith(' '):
                    line_number += 1

            # Get the entire new file.
            new_contents = subprocess.check_output(
                ['git', 'show', new_blob],
                stderr=subprocess.STDOUT).decode("utf-8")

            # Check for other bad things in the entire file.
            if '\r' in new_contents:
                print '      \\r                    {}'.format(new_filename)
                result = 1
            if new_contents.endswith('\n\n'):
                print '      Trailing newlines     {}'.format(new_filename)
                result = 1
            if not new_contents.endswith('\n'):
                print '      Missing last newline  {}'.format(new_filename)
                result = 1

        if ancestor is None:
            # Look for untracked files in important directories
            status_lines = subprocess.check_output(
                ['git', 'status', '-u', 'data', 'horton', 'doc', 'scripts', 'tools',
                 '--porcelain'],
                stderr=subprocess.STDOUT).decode("utf-8").splitlines()
            for status_line in status_lines:
                if status_line.startswith('??'):
                    new_filename = status_line[3:]
                    print '   Untracked file        {}'.format(new_filename)
                    result = 1

    # Stop process with appropriate exit code.
    if result != 0:
        print 'Commit failed. Please, clean up and try again.'
    else:
        print 'OK.'
    sys.exit(result)


if __name__ == '__main__':
    main()
