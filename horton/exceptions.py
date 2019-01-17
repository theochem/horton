# -*- coding: utf-8 -*-
# HORTON: Helpful Open-source Research TOol for N-fermion systems.
# Copyright (C) 2011-2019 The HORTON Development Team
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
'''Definition of all excpetions in HORTON'''


class SymmetryError(Exception):
    '''Exception raised when some symmetry algorithm fails'''
    pass


class ElectronCountError(ValueError):
    '''Exception raised when a negative number of electron is encountered, or
       when more electrons than basis functions are requested.
    '''
    pass


class NoSCFConvergence(Exception):
    '''Exception raised when an SCF algorithm does not reach the convergence
       threshold in the specified number of iterations'''
    pass
