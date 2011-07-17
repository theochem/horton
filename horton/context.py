# -*- coding: utf-8 -*-
# Horton is a Density Functional Theory program.
# Copyright (C) 2011 Toon Verstraelen <Toon.Verstraelen@UGent.be>, ...
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
"""Defines the context in which Horton is used."""


import os
from glob import glob


__all__ = ['context']


class Context(object):
    """Find out where the data directory is located etc.

       The data directory contains data files with standard basis sets and
       pseudo potentials.
    """
    def __init__(self):
        self.data_dir = os.getenv('HORTONDATA')
        if self.data_dir is None:
            self.data_dir = './data'
        if not os.path.isdir(self.data_dir):
            raise IOError('Can not find the data files. The directory %s does not exist.' % self.data_dir)

    def get_fn(self, filename):
        """Return the full path to the given filename"""
        return os.path.join(self.data_dir, filename)

    def glob(self, pattern):
        """Return all files in the data directory that match the given pattern."""
        return glob(self.get_fn(pattern))


context = Context()
