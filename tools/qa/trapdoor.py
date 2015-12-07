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
import cPickle
import json
import os
import shutil
import subprocess
import sys
import time


__all__ = ['TrapdoorProgram']


class TrapdoorProgram(object):
    def __init__(self, name):
        '''Initialize the trapdoor program

           Parameters
           ----------
           name : str
                  The name of the trapdoor program, e.g. ``'cppcheck'``.
        '''
        # Set attributes
        self.name = name
        # Get the QAWORKDIR. Create if it does not exist yet.
        self.qaworkdir = os.getenv('QAWORKDIR', 'qaworkdir')
        if not os.path.isdir(self.qaworkdir):
            os.makedirs(self.qaworkdir)
        self.trapdoor_config_file = os.path.join(self.qaworkdir, 'trapdoor.cfg')

    def main(self):
        args = self.parse_args()
        print r'~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+      ~~~~~~~~~~~~~~~~~'
        print r'  TRAPDOOR %15s:%-7s                      \          _\( )/_' % (
            self.name, args.mode)
        print r'                                                         \          /(o)\ '
        print r'                                                          +~~~~~~~~~~~~~~~~~~~~'
        if args.mode == 'feature':
            self.initialize()
            self.run_tests(args.mode)
        elif args.mode == 'master':
            self.run_tests(args.mode)
        elif args.mode == 'report':
            self.report(args.noisy)
        print

    def parse_args(self):
        '''Parse command-line arguments'''
        parser = argparse.ArgumentParser(prog=os.path.basename(sys.argv[0]))
        parser.add_argument('mode', choices=['feature', 'master', 'report'])
        parser.add_argument('-n', '--noisy', default=False, action='store_true',
            help='Also print output for problems that did not deteriorate.')
        return parser.parse_args()

    def initialize(self):
        shutil.copy('tools/qa/trapdoor.cfg', self.trapdoor_config_file)

    def run_tests(self, mode):
        '''Runs the tests on a single checkout

           Parameters
           ----------
           mode: string
                 A name for the current branch on which the tests are run, typically
                 ``'feature'`` or ``'master'``.

           The results are written to disk in a file ``trapdoor_results_*.pp``. These
           files are later used by the report method to analyze the results.
        '''
        start_time = time.time()
        config = self.load_config()
        counter, messages = self.get_stats(config)
        print 'NUMBER OF MESSAGES', len(messages)
        print 'SUM OF COUNTERS', sum(counter.itervalues())
        fn_pp = 'trapdoor_results_%s_%s.pp' % (self.name, mode)
        with open(os.path.join(self.qaworkdir, fn_pp), 'w') as f:
            cPickle.dump((counter, messages), f)
        print 'WALL TIME %.1f' % (time.time() - start_time)

    def load_config(self):
        with open(self.trapdoor_config_file, 'r') as f:
            return json.load(f)

    def get_stats(self):
        raise NotImplementedError

    def report(self, noisy=False):
        fn_pp_feature = 'trapdoor_results_%s_feature.pp' % self.name
        with open(os.path.join(self.qaworkdir, fn_pp_feature)) as f:
            results_feature = cPickle.load(f)
        fn_pp_master = 'trapdoor_results_%s_master.pp' % self.name
        with open(os.path.join(self.qaworkdir, fn_pp_master)) as f:
            results_master = cPickle.load(f)
        if noisy:
            self.print_details(results_feature, results_master)
        self.check_regression(results_feature, results_master)

    def print_details(self, (counter_feature, messages_feature), (counter_master, messages_master)):
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

        resolved_counter = counter_master - counter_feature
        if len(resolved_counter) > 0:
            print 'SOME COUNTERS DECREASED'
            for key, counter in resolved_counter.iteritems():
                print '%s  |  %+6i' % (key, -counter)


    def check_regression(self, (counter_feature, messages_feature), (counter_master, messages_master)):
        '''Check if the counters got worse

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
        new_messages = sorted(messages_feature - messages_master)
        if len(new_messages) > 0:
            print 'NEW MESSAGES'
            for msg in new_messages:
                print msg

        new_counter = counter_feature - counter_master
        if len(new_counter) > 0:
            print 'SOME COUNTERS INCREASED'
            for key, counter in new_counter.iteritems():
                print '%s  |  %+6i' % (key, counter)
            sys.exit(-1)
        else:
            print 'GOOD ENOUGH'
