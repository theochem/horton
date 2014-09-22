# -*- coding: utf-8 -*-
# Horton is a development platform for electronic structure methods.
# Copyright (C) 2011-2013 Toon Verstraelen <Toon.Verstraelen@UGent.be>
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
'''Two-, three- and four-dimensional matrix implementations

   The purpose of this package is to provide a generic API for different
   implementations of real-valued double precision matrix storage and
   operations.

   Two-dimensional matrices are supposed to be symmetric and are used to
   represent two-index operators and 1DRDMs. Four-dimensional matrices are used
   to represent four-index operators, which are invariant under the following
   interchanges of indexes::

            <ij|kl> = <ji|lk> = <kl|ij> = <lk|ji> =
            <il|kj> = <jk|li> = <kj|il> = <li|jk>

   This module assumes physicists notation for the two-particle operators. It is
   up to the specific implementations of the matrices to make use of these
   symmetries.

   One should use these matrix implementations without accessing the internals
   of each class, i.e. without accessing attributes or methods that start with
   an underscore.

   In order to avoid temporaries when working with arrays, the methods do
   not return arrays. Instead such methods are an in place operation or have
   output arguments. This forces the user to allocate all memory in advance,
   which can then be moved out of the loops. The initial implementation (the
   Dense... classes) are just a proof of concept and may therefore contain
   internals that still make temporaries. This fixed later with an alternative
   implementation.
'''


from horton.matrix.base import *
from horton.matrix.cext import *
from horton.matrix.dense import *
from horton.matrix.cholesky import *
