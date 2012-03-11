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
'''Define a molecular system and all aspects relevant for a computation.

   Objects of the System class specify the geometry (and atomic elements) of the
   molecule on which Horton will perform a computation. Also other parameters
   that determine several aspects of the molecular wavefunction, are attributes
   of a :class:`System` instance, e.g. basissets, pseudopotentials, ghost atoms,
   etc.
'''


import numpy as np

from horton.io import load_system_args


__all__ = ['System']


class System(object):
    def __init__(self, coordinates, numbers, basis=None):
        """
           **Arguments:**

           coordinates
                A (N, 3) float numpy array with Cartesian coordinates of the
                atoms.

           numbers
                A (N,) int numpy vector with the atomic numbers.

           **Optional arguments:**

           basis
                A string or an instance of either the basis set or basis set
                description classes, e.g. 'STO-3G', GOBasisDesc('STO-3G'), ...
        """
        self._coordinates = np.array(coordinates, dtype=float, copy=False)
        self._numbers = np.array(numbers, dtype=int, copy=False)
        from horton.gobasis import GOBasisDesc, GOBasis
        if isinstance(basis, str):
            basis_desc = GOBasisDesc(basis)
            self._basis = basis_desc.apply_to(self)
        elif isinstance(basis, GOBasisDesc):
            self._basis = basis.apply_to(self)
        elif isinstance(basis, GOBasis) or basis is None:
            self._basis = basis
        else:
            raise TypeError('Can not interpret %s as a basis.' % basis)

    def get_natom(self):
        return len(self.numbers)

    natom = property(get_natom)

    def get_coordinates(self):
        return self._coordinates.view()

    coordinates = property(get_coordinates)

    def get_numbers(self):
        return self._numbers.view()

    numbers = property(get_numbers)

    def get_basis(self):
        return self._basis

    basis = property(get_basis)

    @classmethod
    def from_file(cls, *args, **kwargs):
        """Create a System object from a file.

           A list of filenames may be provided, which will be loaded in that
           order. Each file complements or overrides the information loaded
           from a previous file in the list. Furthermore, keyword arguments
           may be used to specify additional constructor arguments.
        """
        constructor_args = {}
        for fn in args:
            fn_args = load_system_args(fn)
            constructor_args.update(fn_args)
        constructor_args.update(kwargs)
        return cls(**constructor_args)
