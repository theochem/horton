# -*- coding: utf-8 -*-
# HORTON: Helpful Open-source Research TOol for N-fermion systems.
# Copyright (C) 2011-2016 The HORTON Development Team
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
import bisect
import cPickle
from fnmatch import fnmatch
import json
import os
import shutil
import subprocess
import sys
import time


__all__ = ['Message', 'TrapdoorProgram', 'get_source_filenames', 'run_command']


class Message(object):
    """Error message and meta information.

    This class contains all the machinery that the trapdoor driver uses to decide if
    a message is new or not. When a context is available, it is used to compare two
    messages instead of line numbers. If not, the line numbers are used. Line numbers are
    a relatively poor descriptor for assessing if a message is new in a feature branch.
    For example, lines may have been inserted or removed, which changes the line numbers
    without actually changing any code.
    """

    def __init__(self, filename, lineno, charno, text, context=None):
        """Initialize a message.

        Parameters
        ----------
        filename : str
                   The filename for which the message is reported.
        lineno : int (or None)
                 The line number at which the error is reported, if any.
        charno : int (or None)
                 The character position at which the error is reported, if any.
        text : str
               A description of the problem.
        context : str
                  A string that (almost) uniquely identifies the location of the error
                  without using line numbers.
        """
        if lineno is not None and not isinstance(lineno, int):
            raise TypeError('When given, lineno must be integer.')
        if charno is not None and not isinstance(charno, int):
            raise TypeError('When given, charno must be integer.')
        self._filename = filename
        self._lineno = lineno
        self._charno = charno
        self._text = text
        self._context = context

    filename = property(lambda self: self._filename)
    lineno = property(lambda self: self._lineno)
    charno = property(lambda self: self._charno)
    text = property(lambda self: self._text)
    context = property(lambda self: self._context)

    def __eq__(self, other):
        """Test if self is equal to other."""
        # First come the usualy things for comparisons...
        equal = self.__class__ == other.__class__ \
            and self._filename == other._filename \
            and self._charno == other._charno \
            and self._text == other._text \
            and self._context == other._context
        if not equal:
            return False
        # If still equal, then only use line numbers if no context is available.
        if self._context is None and other._context is None:
            return self._lineno == other._lineno
        return True

    def __hash__(self):
        """Return a fast hash.

        The hash only includes the lineno if no context is present. If a context is
        present, it is used instead and the line number is not included in the hash. This
        convention is compatible with the code in __eq__.
        """
        if self._context is None:
            return hash((self._filename, self._lineno, self._charno, self._text))
        else:
            return hash((self._filename, self._charno, self._text, self._context))

    def __lt__(self, other):
        """Test if self is less than other."""
        if self.__class__ != other.__class__:
            return self < other
        tup_self = (self._filename, self._lineno, self._charno, self._text, self._context)
        tup_other = (other._filename, other._lineno, other._charno, other._text, other._context)
        return tup_self < tup_other

    def add_context(self, context):
        """Return an identical message with context."""
        if self._context is not None:
            raise ValueError('This message already has context.')
        return Message(self._filename, self._lineno, self._charno, self._text, context)

    def __str__(self):
        """Return a nicely formatted string representation of the message."""
        # Fix the location string
        if self.filename is None:
            location = '(nofile)            '
        else:
            location = str(self.filename)
            if self.lineno is not None:
                location += '%6i' % self.lineno
            else:
                location += ' '*6
            if self.charno is not None:
                location += '%6i' % self.charno
            else:
                location += ' '*6
        return '%70s   %s' % (location, self.text)


def _print_messages(header, messages, pattern=None):
    """Print a set of messages.

    Parameters
    ----------
    header : str
             A Header string, usually uppercase.
    messages : iterable
               A list of Message instances.
    pattern : None or str
              When given, only messages containing ``pattern`` will be printed.
    """
    if len(messages) > 0:
        print header
        for msg in sorted(messages):
            if pattern is None or pattern in msg.filename:
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
                            help='Only print messages whose filename contains PATTERN')
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
        print 'ADDING SOURCE ...'
        self._add_contexts(messages)
        print 'NUMBER OF MESSAGES :', len(messages)
        print 'SUM OF COUNTERS    :', sum(counter.itervalues())
        fn_pp = 'trapdoor_results_%s_%s.pp' % (self.name, mode)
        with open(os.path.join(self.qaworkdir, fn_pp), 'w') as f:
            cPickle.dump((counter, messages), f)
        print 'WALL TIME          : %.1f' % (time.time() - start_time)

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
        messages : Set([]) of Message instances
                   All errors encountered in the current branch.
        """
        raise NotImplementedError

    def _add_contexts(self, all_messages):
        """Add source lines to the messages.

        This method has the desired side effect that messages that only differ in line
        number but that do have the same source line, will be considered identical.

        Parameters
        ----------
        all_messages : Set([]) of Message instances
                       All errors encountered in the current branch.
        """
        # 1) Collect all messages in a dictionary, where filenames are keys and values
        #    are sorted lists of messages. Only messages with a filename and a line number
        #    must be included.
        mdict = {}
        for message in all_messages:
            if message.filename is not None and message.lineno is not None:
                l = mdict.setdefault(message.filename, [])
                bisect.insort(l, message)
        # 2) Loop over all files and collect some source context for each message
        for filename, file_messages in mdict.iteritems():
            with open(filename) as source_file:
                lines = source_file.readlines()
                for message in file_messages:
                    all_messages.discard(message)
                    # The context starts three lines before the line and ends three lines
                    # after.
                    context = ''.join(lines[
                        max(0, message.lineno - 3):
                        min(len(lines), message.lineno + 4)
                    ])
                    all_messages.add(message.add_context(context))

    def report(self, noisy=False, pattern=None):
        """Load feature and ancestor results from disk and report on screen.

        Parameters
        ----------
        noisy : bool
                If True, more detailed screen output is printed.
        pattern : None or str
                  When given, only messages containing ``pattern`` will be printed.
        """
        # Load all the trapdoor results from the two branches.
        fn_pp_feature = 'trapdoor_results_%s_feature.pp' % self.name
        with open(os.path.join(self.qaworkdir, fn_pp_feature)) as f:
            results_feature = cPickle.load(f)
        fn_pp_ancestor = 'trapdoor_results_%s_ancestor.pp' % self.name
        with open(os.path.join(self.qaworkdir, fn_pp_ancestor)) as f:
            results_ancestor = cPickle.load(f)
        # Make the report.
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
        messages_feature : Set([]) of Message instances
                           All errors encountered in the feature branch.
        counter_ancestor : collections.Counter
                           Counts for different error types in the ancestor.
        messages_ancestor : Set([]) of Message instances
                            All errors encountered in the ancestor.
        pattern : None or str
                  When given, only messages containing ``pattern`` will be printed.
        """
        resolved_messages = sorted(messages_ancestor - messages_feature)
        _print_messages('RESOLVED MESSAGES', resolved_messages, pattern)

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
        messages_feature : Set([]) of Message instances
                           All errors encountered in the feature branch.
        counter_ancestor : collections.Counter
                         Counts for different error types in the ancestor.
        messages_ancestor : Set([]) of Message instances
                          All errors encountered in the ancestor.
        pattern : None or str
                  When given, only messages containing ``pattern`` will be printed.
        """
        new_messages = sorted(messages_feature - messages_ancestor)
        _print_messages('NEW MESSAGES', new_messages, pattern)

        new_counter = counter_feature - counter_ancestor
        if len(new_counter) > 0:
            print 'SOME COUNTERS INCREASED'
            for key, counter in new_counter.iteritems():
                print '%s  |  %+6i' % (key, counter)
            sys.exit(1)
        else:
            print 'GOOD (ENOUGH)'


def get_source_filenames(config, language, unpackaged_only=False):
    """Return a list of source files according to configuration settings.

    This function will search for all cpp or py source files in the "XXX_directories" and
    will exclude some based in fnmatch patterns provided in the "XXX_exclude" setting,
    where XXX is 'cpp' or 'py'.

    Parameters
    ----------
    config : dict
             A dictionary of configuration settings, loaded with json from trapdoor.cfg.
    language : str
               'py' or 'cpp'
    unpackaged_only : bool
                      When set to True, only files are listed that are not part of a
                      (Python) package.
    """
    # Extensions for each language:
    if language == 'cpp':
        source_patterns = ['*.cpp', '*.h']
    elif language == 'py':
        source_patterns = ['*.py', '*.pyx']
    else:
        raise ValueError('Language must be either \'cpp\' or \'py\'.')

    # Get config settings
    directories = config['%s_directories' % language]
    exclude = config['%s_exclude' % language]
    packages = config.get('%s_packages' % language, [])

    # Define the filename filter
    def acceptable(dirpath, filename):
        """Determine if a filename should be included in the result."""
        if not any(fnmatch(filename, source_pattern) for source_pattern in source_patterns):
            return False
        if not all(not fnmatch(filename, exclude_filter) for exclude_filter in exclude):
            return False
        if unpackaged_only:
            in_package = False
            for package in packages:
                if dirpath.startswith(package):
                    in_package = True
                    break
            if in_package:
                return False
        return True

    # Loop over all files in given directories
    result = []
    for source_directory in directories:
        for dirpath, _dirnames, filenames in os.walk(source_directory):
            for filename in filenames:
                if acceptable(dirpath, filename):
                    result.append(os.path.join(dirpath, filename))
    return result


def run_command(command, verbose=True, cwd=None, has_failed=None):
    """Run command as subprocess with default settings suitable for trapdoor scripts.

    Parameters
    ----------
    command : list of str
              The command argument to be passed to Popen.
    verbose : bool
              When set to False, the command will not be printed on screen.
    cwd : str
          The working directory where the command is executed.
    has_failed : function(returncode, stdout, stderr)
                 A function that determines if the subprocess has failed. The default
                 behavior is to check for a non-zero return code.

    Returns
    -------
    output : (str, str)
             if Sucessful, the output collected from stdout and stderr are returned.

    Raises
    ------
    In case the subprocess returns a non-zero exit code, the stdout and stderr are printed
    on screen and RuntimeError is raised.
    """
    # Functions to detect failure
    def default_has_failed(returncode, stdout, stderr):
        """Default function to detect failed subprocess."""
        return returncode != 0
    if has_failed is None:
        has_failed = default_has_failed

    if verbose:
        print 'RUNNING            :', ' '.join(command)
    proc = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=cwd)
    stdout, stderr = proc.communicate()
    if has_failed(proc.returncode, stdout, stderr):
        print 'STDOUT'
        print '------'
        print stdout
        print 'STDERR'
        print '------'
        print stderr
        raise RuntimeError('Subprocess returned non-zero exit status %i' % proc.returncode)
    else:
        return stdout, stderr
