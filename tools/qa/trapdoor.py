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
"""Shared trapdoor code.

This model provides the ``TrapdoorProgram`` base class for all trapdoor programs.
"""


import argparse
import cPickle
import json
import os
import shutil
import sys
import time


__all__ = ['TrapdoorProgram']


def _print_messages(header, messages, pattern=None):
    """Print a set of messages.

    Parameters
    ----------
    header : str
             A Header string, usually uppercase
    messages : iterable
               A list of message strings
    pattern : None or str
              When given, only messages containing ``pattern`` will be printed.
    """
    if len(messages) > 0:
        print header
        for msg in messages:
            if pattern is None or pattern in msg:
                print msg


class TrapdoorProgram(object):
    """Base class for all trapdoor programs.

    This class implements all the shared functionality between different specializations
    of the traps door programs, such as the command-line interface, basic configuration
    and some screen output.

    Trapdoor programs must be implemented as follows:

    * Extend the ``__init__`` method, at least providing a sensible name.
    * Optional extend the ``prepare`` method, e.g. to copy config files.
    * Override the ``get_stats`` method, where the real work is done: calling a QA program
      and collecting its output into a standard format.
    * Call ``DerivedTrapdoorProgram().main()``
    """

    def __init__(self, name):
        """Initialize the trapdoor program.

        Parameters
        ----------
        name : str
               The name of the trapdoor program, e.g. ``'cppcheck'``.
        """
        # Set attributes
        self.name = name
        # Get the QAWORKDIR. Create if it does not exist yet.
        self.qaworkdir = os.getenv('QAWORKDIR', 'qaworkdir')
        if not os.path.isdir(self.qaworkdir):
            os.makedirs(self.qaworkdir)
        self.trapdoor_config_file = os.path.join(self.qaworkdir, 'trapdoor.cfg')

    def main(self):
        """Execute the main routine of the trapdoor program.

        This includes parsing command-line arguments and running one of the three modes:
        ``feature``, ``ancestor`` or ``report``.
        """
        args = self.parse_args()
        print r'~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+      ~~~~~~~~~~~~~~~~~'
        print r'  TRAPDOOR %15s:%-10s                   \          _\( )/_' % (
            self.name, args.mode)
        print r'                                                         \          /(o)\ '
        print r'                                                          +~~~~~~~~~~~~~~~~~~~~'
        if args.mode == 'feature':
            self.prepare()
            self.run_tests(args.mode)
        elif args.mode == 'ancestor':
            self.run_tests(args.mode)
        elif args.mode == 'report':
            self.report(args.noisy, args.pattern)
        print

    def parse_args(self):
        """Parse command-line arguments.

        Returns
        -------
        args : argsparse.Namespace
               The parsed command-line arguments.
        """
        parser = argparse.ArgumentParser(prog=os.path.basename(sys.argv[0]))
        parser.add_argument('mode', choices=['feature', 'ancestor', 'report'])
        parser.add_argument('-n', '--noisy', default=False, action='store_true',
                            help='Also print output for problems that did not '
                                 'deteriorate.')
        parser.add_argument('-f', '--pattern', metavar='PATTERN', dest='pattern',
                            help='Only print messages containing PATTERN')
        return parser.parse_args()

    def prepare(self):
        """Prepare for the tests, only once, needed for both feature branch and ancestor.

        This usually comes down to copying some config files to the QAWORKDIR. This method
        is only called when in the feature branch.
        """
        shutil.copy('tools/qa/trapdoor.cfg', self.trapdoor_config_file)

    def run_tests(self, mode):
        """Run the tests on a single branch HEAD.

        Parameters
        ----------
        mode: string
              A name for the current branch on which the tests are run, typically
              ``'feature'`` or ``'ancestor'``.

        The results are written to disk in a file ``trapdoor_results_*.pp``. These files
        are later used by the report method to analyze the results.
        """
        start_time = time.time()
        with open(self.trapdoor_config_file, 'r') as f:
            config = json.load(f)
        counter, messages = self.get_stats(config)
        print 'NUMBER OF MESSAGES :', len(messages)
        print 'SUM OF COUNTERS    :', sum(counter.itervalues())
        fn_pp = 'trapdoor_results_%s_%s.pp' % (self.name, mode)
        with open(os.path.join(self.qaworkdir, fn_pp), 'w') as f:
            cPickle.dump((counter, messages), f)
        print 'WALL TIME %.1f' % (time.time() - start_time)

    def get_stats(self, config):
        """Run tests using an external program and collect its output.

        This method must be implemented in a subclass.

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
        raise NotImplementedError

    def report(self, noisy=False, pattern=None):
        """Load feature and ancestor results from disk and report on screen.

        Parameters
        ----------
        noisy : bool
                If True, more detailed screen output is printed.
        pattern : None or str
                  When given, only messages containing ``pattern`` will be printed.
        """
        fn_pp_feature = 'trapdoor_results_%s_feature.pp' % self.name
        with open(os.path.join(self.qaworkdir, fn_pp_feature)) as f:
            results_feature = cPickle.load(f)
        fn_pp_ancestor = 'trapdoor_results_%s_ancestor.pp' % self.name
        with open(os.path.join(self.qaworkdir, fn_pp_ancestor)) as f:
            results_ancestor = cPickle.load(f)
        if noisy:
            self.print_details(results_feature, results_ancestor, pattern)
        self.check_regression(results_feature, results_ancestor, pattern)

    def print_details(self, (counter_feature, messages_feature),
                      (counter_ancestor, messages_ancestor), pattern=None):
        """Print optional detailed report of the test results.

        Parameters
        ----------
        counter_feature : collections.Counter
                          Counts for different error types in the feature branch.
        messages_feature : Set([]) of strings
                           All errors encountered in the feature branch.
        counter_ancestor : collections.Counter
                         Counts for different error types in the ancestor.
        messages_ancestor : Set([]) of strings
                          All errors encountered in the ancestor.
        pattern : None or str
                  When given, only messages containing ``pattern`` will be printed.
        """
        resolved_messages = sorted(messages_ancestor - messages_feature)
        _print_messages('RESOLVED MESSAGES (also includes messages for which just the '
                        'line number changed)', resolved_messages, pattern)

        unchanged_messages = sorted(messages_ancestor & messages_feature)
        _print_messages('UNCHANGED MESSAGES', unchanged_messages, pattern)

        resolved_counter = counter_ancestor - counter_feature
        if len(resolved_counter) > 0:
            print 'SOME COUNTERS DECREASED'
            for key, counter in resolved_counter.iteritems():
                print '%s  |  %+6i' % (key, -counter)

    def check_regression(self, (counter_feature, messages_feature),
                         (counter_ancestor, messages_ancestor), pattern=None):
        """Check if the counters got worse.

        The new errors are printed and if a regression is observed, the program quits with
        an exit status of 1.

        Parameters
        ----------
        counter_feature : collections.Counter
                          Counts for different error types in the feature branch.
        messages_feature : Set([]) of strings
                           All errors encountered in the feature branch.
        counter_ancestor : collections.Counter
                         Counts for different error types in the ancestor.
        messages_ancestor : Set([]) of strings
                          All errors encountered in the ancestor.
        pattern : None or str
                  When given, only messages containing ``pattern`` will be printed.
        """
        new_messages = sorted(messages_feature - messages_ancestor)
        _print_messages('NEW MESSAGES (also includes messages for which just the line '
                        'number changed)', new_messages, pattern)

        new_counter = counter_feature - counter_ancestor
        if len(new_counter) > 0:
            print 'SOME COUNTERS INCREASED'
            for key, counter in new_counter.iteritems():
                print '%s  |  %+6i' % (key, counter)
            sys.exit(1)
        else:
            print 'GOOD (ENOUGH)'
