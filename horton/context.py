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
'''The context in which HORTON is used

   This module controls global parameters that are purely technical, e.g. the
   location of data files. It is certainly not meant to keep track of input
   parameters for a computation.

   This module contains a context object, an instance of the :class:`Context`
   class. For now, its functionality is rather limited. It tries to figure
   out the location of the data directory. If it is not specified in the
   environment variable ``HORTONDATA``, it is assumed that the data is located
   in a directory called ``data``. If the data directory does not exist, an
   error is raised.
'''


import os, sys
from glob import glob


__all__ = ['context', 'Context']


class Context(object):
    '''Finds out where the data directory is located etc.

       The data directory contains data files with standard basis sets and
       pseudo potentials.
    '''
    def __init__(self):
        # Determine data directory (also for in-place build)
        self.data_dir = os.getenv('HORTONDATA')
        if self.data_dir is None:
            fn_data_dir = os.path.join(os.path.dirname(__file__), 'data_dir.txt')
            if os.path.isfile(fn_data_dir):
                with open(fn_data_dir) as f:
                    self.data_dir = os.path.join(f.read().strip(), 'share/horton')
        if self.data_dir is None:
            self.data_dir = './data'
        self.data_dir = os.path.abspath(self.data_dir)
        # Determine include directory
        self.include_dir = os.getenv('HORTONINCLUDE')
        if self.include_dir is None:
            fn_data_dir = os.path.join(os.path.dirname(__file__), 'data_dir.txt')
            if os.path.isfile(fn_data_dir):
                with open(fn_data_dir) as f:
                    self.include_dir = os.path.join(
                        f.read().strip(),
                        'include/python%i.%i' % (sys.version_info.major, sys.version_info.minor))
        if not os.path.isdir(self.data_dir):
            raise IOError('Can not find the data files. The directory %s does not exist.' % self.data_dir)

    def get_fn(self, filename):
        '''Return the full path to the given filename in the data directory.'''
        return os.path.join(self.data_dir, filename)

    def glob(self, pattern):
        '''Return all files in the data directory that match the given pattern.'''
        return glob(self.get_fn(pattern))

    def get_include(self):
        '''Return the list with directories containing header files (.h and .pxd)'''
        return self.include_dir


context = Context()
