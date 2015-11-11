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
'''Trapdoor driver routines

   The ``main`` function in this script runs tests on the feature and on the master
   branch. It can be used by other trapdoor_*.py scripts as main routine.
'''


import argparse
import os
import subprocess
import sys


__all__ = ['main']


def main(get_stats):
    '''Main routine doing all the work for a given code quality check

       Parameters
       ----------
       get_stats : function
                   This function takes no parameters and returns ``counter`` and
                   ``messages``. It does quality checks on the source code that is
                   currently checked out. See ``trapdoor_cppcheck.py`` for an example. The
                   ``counter`` return value must be of the type ``collections.Counter``.
                   It counts all types of error messages. The idea is that these counts
                   should not be allowed to increase. Please use relatively specific keys
                   in ``counter``, e.g. containing filename and the type of error. The
                   ``messages`` return value is a ``Set`` of strings with all the
                   encountered error messages.

    '''
    args = parse_args()
    counter_feature, messages_feature, counter_master, messages_master = run_tests(get_stats)
    if args.noisy:
        print_details(counter_feature, messages_feature, counter_master, messages_master)
    check_deterioration(counter_feature, counter_master)


def parse_args():
    '''Parse command-line arguments'''
    parser = argparse.ArgumentParser(prog=os.path.basename(sys.argv[0]))
    parser.add_argument('-n', '--noisy', default=False, action='store_true',
        help='Also print output for problems that did not deteriorate.')
    return parser.parse_args()


def run_tests(get_stats):
    '''Runs the tests on two checkouts and returns counters and messages

       Parameters
       ----------
       get_stats: function, no parameters, returns ``counter`` and ``messages``
                  A function that does quality checks on the source code that is currently
                  checked out. See ``trapdoor_cppcheck.py`` for an example.

       Returns
       -------
       counter_feature: collections.Counter
                        counts for different error types in the feature branch
       messages_feature: Set([]) of strings
                         all errors encountered in the feature branch
       counter_master: collections.Counter
                        counts for different error types in the master branch
       messages_master: Set([]) of strings
                         all errors encountered in the master branch
    '''
    # Get git info
    name_feature = subprocess.check_output(['git', 'rev-parse', '--abbrev-ref', 'HEAD']).strip()
    commit_id_feature = subprocess.check_output(['git', 'rev-parse', 'HEAD']).strip()
    commit_id_master = subprocess.check_output(['git', 'rev-parse', 'master']).strip()
    # Actual work
    counter_feature, messages_feature = get_stats()
    print 'CHECKING OUT master (%s)' % (commit_id_master)
    subprocess.call(['git', 'checkout', 'master'])
    counter_master, messages_master = get_stats()
    print 'CHECKING OUT %s (%s)' % (name_feature, commit_id_feature)
    subprocess.call(['git', 'checkout', name_feature])
    return counter_feature, messages_feature, counter_master, messages_master


def print_details(counter_feature, messages_feature, counter_master, messages_master):
    '''Print optional detailed report of the test results

       Parameters
       ----------
       counter_feature: collections.Counter
                        counts for different error types in the feature branch
       messages_feature: Set([]) of strings
                         all errors encountered in the feature branch
       counter_master: collections.Counter
                        counts for different error types in the master branch
       messages_master: Set([]) of strings
                         all errors encountered in the master branch
    '''
    resolved_messages = sorted(messages_master - messages_feature)
    if len(resolved_messages) > 0:
        print 'RESOLVED MESSAGES'
        for msg in resolved_messages:
            print msg

    unchanged_messages = sorted(messages_master & messages_feature)
    if len(unchanged_messages) > 0:
        print 'UNCHANGED MESSAGES'
        for msg in unchanged_messages:
            print msg

    new_messages = sorted(messages_feature - messages_master)
    if len(new_messages) > 0:
        print 'NEW MESSAGES'
        for msg in new_messages:
            print msg

    resolved_counter = counter_master - counter_feature
    if len(resolved_counter) > 0:
        print 'SOME STATISTICS IMPROVED'
        for key, counter in resolved_counter.iteritems():
            print '%s  |  %+6i' % (key, -counter)


def check_deterioration(counter_feature, counter_master):
    '''Check if the counters got worse

       Parameters
       ----------
       counter_feature: collections.Counter
                        counts for different error types in the feature branch
       counter_master: collections.Counter
                        counts for different error types in the master branch
    '''
    new_counter = counter_feature - counter_master
    if len(new_counter) > 0:
        print 'SOME STATISTICS GOT WORSE'
        for key, counter in new_counter.iteritems():
            print '%s  |  %+6i' % (key, counter)
        sys.exit(-1)
