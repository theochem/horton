# -*- coding: utf-8 -*-
# Horton is a Density Functional Theory program.
# Copyright (C) 2011-2012 Toon Verstraelen <Toon.Verstraelen@UGent.be>
#
# This file is part of Horton.
#
# Horton is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
#
# Horton is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, see <http://www.gnu.org/licenses/>
#
#--
'''Screen logging, timing and citation management.

   The goal of the screen logger is to track the progress of a computation in
   a convenient human-readable way, possibly higlighting problematic situations.
   It is not intended as a computer-readable output file that contains all the
   results of a computation. For that purpose, all useful information is
   written to a binary checkpoint file or kept in memory as attributes of the
   Horton objects.
'''

# TODO: - track memory usage
# TODO: - track references that should be cited
# TODO: - use timer
# TODO: - use the logger in all current code
# TODO: - find more convenient ways of using the logger

import sys, os, datetime, getpass, time, codecs, locale, functools, atexit
from contextlib import contextmanager
import horton


__all__ = ['log', 'timer']


class ScreenLog(object):
    # log levels
    silent = 0
    warning = 1
    low = 2
    medium = 3
    high = 4
    debug = 5

    # screen parameters
    margin = 9
    width = 71

    def __init__(self, name, version, head_banner, foot_banner, timer, f=None):
        self.name = name
        self.version = version
        self.head_banner = head_banner
        self.foot_banner = foot_banner
        self.timer = timer

        self._active = False
        self._level = self.medium
        self.prefix = ' '*(self.margin-1)
        self._last_used_prefix = None
        self.stack = []
        self.add_newline = False
        if f is None:
            _file = sys.stdout
        else:
            _file = f
        # Wrap sys.stdout into a StreamWriter to allow writing unicode.
        self._file = codecs.getwriter(locale.getpreferredencoding())(_file)

    do_warning = property(lambda self: self._level >= self.warning)
    do_low = property(lambda self: self._level >= self.low)
    do_medium = property(lambda self: self._level >= self.medium)
    do_high = property(lambda self: self._level >= self.high)
    do_debug = property(lambda self: self._level >= self.debug)

    def _pass_color_code(self, code):
        if self._file.isatty():
            return code
        else:
            return ''

    reset =     property(functools.partial(_pass_color_code, code="\033[0m"))
    bold =      property(functools.partial(_pass_color_code, code="\033[01m"))
    teal =      property(functools.partial(_pass_color_code, code="\033[36;06m"))
    turquoise = property(functools.partial(_pass_color_code, code="\033[36;01m"))
    fuchsia =   property(functools.partial(_pass_color_code, code="\033[35;01m"))
    purple =    property(functools.partial(_pass_color_code, code="\033[35;06m"))
    blue =      property(functools.partial(_pass_color_code, code="\033[34;01m"))
    darkblue =  property(functools.partial(_pass_color_code, code="\033[34;06m"))
    green =     property(functools.partial(_pass_color_code, code="\033[32;01m"))
    darkgreen = property(functools.partial(_pass_color_code, code="\033[32;06m"))
    yellow =    property(functools.partial(_pass_color_code, code="\033[33;01m"))
    brown =     property(functools.partial(_pass_color_code, code="\033[33;06m"))
    red =       property(functools.partial(_pass_color_code, code="\033[31;01m"))


    def set_level(self, level):
        if level < self.silent or level > self.debug:
            raise ValueError('The level must be one of the ScreenLog attributes.')
        self._level = level

    def __call__(self, *words):
        s = u' '.join(unicode(w) for w in words)
        if not self.do_warning:
            raise RuntimeError('The runlevel should be at least warning when logging.')
        if not self._active:
            prefix = self.prefix
            self.print_header()
            self.prefix = prefix
        if self.add_newline and self.prefix != self._last_used_prefix:
            print >> self._file
            self.add_newline = False
        # Check for alignment code '&'
        pos = s.find(u'&')
        if pos == -1:
            lead = u''
            rest = s
        else:
            lead = s[:pos] + ' '
            rest = s[pos+1:]
        width = self.width - len(lead)
        if width < self.width/2:
            raise ValueError('The lead may not exceed half the width of the terminal.')
        # break and print the line
        first = True
        while len(rest) > 0:
            if len(rest) > width:
                pos = rest.rfind(' ', 0, width)
                if pos == -1:
                    current = rest[:width]
                    rest = rest[width:]
                else:
                    current = rest[:pos]
                    rest = rest[pos:].lstrip()
            else:
                current = rest
                rest = u''
            print >> self._file, u'%s %s%s' % (self.prefix, lead, current)
            if first:
                lead = u' '*len(lead)
                first = False
        self._last_used_prefix = self.prefix

    def warn(self, *words):
        text = u'WARNING!!&'+u' '.join(unicode(w) for w in words)
        text = '%s%s%s' % (self.red, text, self.reset)
        self(text)

    def hline(self, char='~'):
        self(char*self.width)

    def center(self, *words, **kwargs):
        if len(kwargs) == 0:
            edge = ''
        elif len(kwargs) == 1:
            if 'edge' not in kwargs:
                raise TypeError('Only one keyword argument is allowed, that is edge')
            edge = kwargs['edge']
        else:
            raise TypeError('Too many keyword arguments. Should be at most one.')
        s = u' '.join(unicode(w) for w in words)
        if len(s) + 2*len(edge) > self.width:
            raise ValueError('Line too long. center method does not support wrapping.')
        self('%s%s%s' % (edge, s.center(self.width-2*len(edge)), edge))

    def blank(self):
        print >> self._file

    def _enter(self, prefix):
        if len(prefix) > self.margin-1:
            raise ValueError('The prefix must be at most %s characters wide.' % (self.margin-1))
        self.stack.append(self.prefix)
        self.prefix = prefix.upper().rjust(self.margin-1, ' ')
        self.add_newline = True

    def _exit(self):
        self.prefix = self.stack.pop(-1)
        if self._active:
            self.add_newline = True

    @contextmanager
    def section(self, prefix):
        self._enter(prefix)
        try:
            yield
        finally:
            self._exit()

    def print_header(self):
        if self.do_warning and not self._active:
            self._active = True
            print >> self._file, self.head_banner
            self._print_basic_info()

    def print_footer(self):
        if self.do_warning and self._active:
            self._print_basic_info()
            self.timer._stop('Total')
            self.timer.report(self)
            print >> self._file, self.foot_banner

    def _print_basic_info(self):
        if self.do_low:
            with self.section('ENV'):
                self('User:          &' + getpass.getuser())
                self('Machine info:  &' + ' '.join(os.uname()))
                self('Time:          &' + datetime.datetime.now().isoformat())
                self('Python version:&' + sys.version.replace('\n', ''))
                self('%s&%s' % (('%s version:' % self.name).ljust(15), self.version))
                self('Current Dir:   &' + os.getcwd())
                self('Command line:  &' + ' '.join(sys.argv))


class Timer(object):
    def __init__(self):
        self.cpu = 0.0
        self._start = None

    def start(self):
        assert self._start is None
        self._start = time.clock()

    def stop(self):
        assert self._start is not None
        self.cpu += time.clock() - self._start
        self._start = None


class SubTimer(object):
    def __init__(self, label):
        self.label = label
        self.total = Timer()
        self.own = Timer()

    def start(self):
        self.total.start()
        self.own.start()

    def start_sub(self):
        self.own.stop()

    def stop_sub(self):
        self.own.start()

    def stop(self):
        self.own.stop()
        self.total.stop()


class TimerGroup(object):
    def __init__(self):
        self.parts = {}
        self._stack = []
        self._start('Total')

    def reset(self):
        for timer in self.parts.itervalues():
            timer.total.cpu = 0.0
            timer.own.cpu = 0.0

    @contextmanager
    def section(self, label):
        self._start(label)
        try:
            yield
        finally:
            self._stop(label)

    def _start(self, label):
        # get the right timer object
        timer = self.parts.get(label)
        if timer is None:
            timer = SubTimer(label)
            self.parts[label] = timer
        # start timing
        timer.start()
        if len(self._stack) > 0:
            self._stack[-1].start_sub()
        # put it on the stack
        self._stack.append(timer)

    def _stop(self, label):
        timer = self._stack.pop(-1)
        assert timer.label == label
        timer.stop()
        if len(self._stack) > 0:
            self._stack[-1].stop_sub()

    def get_max_own_cpu(self):
        result = None
        for part in self.parts.itervalues():
            if result is None or result < part.own.cpu:
                result = part.own.cpu
        return result

    def report(self, log):
        max_own_cpu = self.get_max_own_cpu()
        #if max_own_cpu == 0.0:
        #    return
        with log.section('TIMER'):
            log('Overview of CPU time usage.')
            log.hline()
            log('Label             Total      Own')
            log.hline()
            bar_width = log.width-33
            for label, timer in sorted(self.parts.iteritems()):
                #if timer.total.cpu == 0.0:
                #    continue
                if max_own_cpu > 0:
                    cpu_bar = "W"*int(timer.own.cpu/max_own_cpu*bar_width)
                else:
                    cpu_bar = ""
                log('%14s %8.1f %8.1f %s' % (
                    label.ljust(14),
                    timer.total.cpu, timer.own.cpu, cpu_bar.ljust(bar_width),
                ))
            log.hline()


head_banner = """\
 _ __ _
/ |..| \ Welcome to Horton, an open source constrained HF/DFT/1RDM program.
\/ || \/
 |_''_|  Horton is written by Toon Verstraelen (1) and Matthew Chan (2).

         (1) Center for Molecular Modeling, Ghent University, Belgium.
         (2) The Ayers Group at McMaster University, Ontario, Canada.

         More information about Horton can be found on this website:
         http://theochem.github.com/horton/

         The purpose of this log file is to track the progress and quality of a
         computation. Useful numerical output may be written to a checkpoint
         file and is accessible through the Python scripting interface.

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"""


foot_banner = """
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 _    _
/ >--< \ End of the Horton program.
\|  \ |/
 |_||_|  Thank you for using Horton! See you soon!
"""

timer = TimerGroup()
log = ScreenLog('HORTON', horton.__version__, head_banner, foot_banner, timer)
atexit.register(log.print_footer)
