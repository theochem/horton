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
from horton.matrix import DenseLinalgFactory


__all__ = ['System']


class System(object):
    def __init__(self, coordinates, numbers, basis=None, wfn=None, lf=None, operators=None):
        """
           **Arguments:**

           coordinates
                A (N, 3) float numpy array with Cartesian coordinates of the
                atoms.

           numbers
                A (N,) int numpy vector with the atomic numbers.

           basis
                A string or an instance of either the basis set or basis set
                description classes, e.g. 'STO-3G', GOBasisDesc('STO-3G'), ...

           **Optional arguments:**

           wfn
                A wavefunction object.

           lf
                A LinalgFactory instance. When not given, a DenseLinalgFactory
                is used by default.

           operators
                A dictionary with one- and two-body operators.
        """
        # A) Assign all attributes
        self._coordinates = np.array(coordinates, dtype=float, copy=False)
        self._numbers = np.array(numbers, dtype=int, copy=False)
        #
        from horton.gbasis import GOBasisDesc, GOBasis
        if isinstance(basis, str):
            basis_desc = GOBasisDesc(basis)
            self._basis = basis_desc.apply_to(self)
        elif isinstance(basis, GOBasisDesc):
            self._basis = basis.apply_to(self)
        elif isinstance(basis, GOBasis) or basis is None:
            self._basis = basis
        else:
            raise TypeError('Can not interpret %s as a basis.' % basis)
        #
        self._wfn = wfn
        if wfn is not None:
            assert wfn.basis == basis
        #
        if lf is None:
            self._lf = DenseLinalgFactory()
        else:
            self._lf = lf
        #
        if operators is None:
            self._operators = {}
        else:
            self._operators = operators

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

    def get_wfn(self):
        return self._wfn

    wfn = property(get_wfn)

    def get_lf(self):
        return self._lf

    lf = property(get_lf)

    def get_operators(self):
        return self._operators

    operators = property(get_operators)

    @classmethod
    def from_file(cls, *args, **kwargs):
        """Create a System object from a file.

           A list of filenames may be provided, which will be loaded in that
           order. Each file complements or overrides the information loaded
           from a previous file in the list. Furthermore, keyword arguments
           may be used to specify additional constructor arguments.

           The ``lf`` optional argument is picked up from the kwargs list to
           contstruct (when needed) arrays to store the results loaded from
           file. When ``lf`` is not given, a DenseLinalgFactory is created by
           default.
        """
        constructor_args = {}
        lf = kwargs.get('lf')
        if lf is None:
            lf = DenseLinalgFactory()
        for fn in args:
            fn_args = load_system_args(fn, lf)
            constructor_args.update(fn_args)
        constructor_args.update(kwargs)

        # If the basis comes from an external code and some operators are
        # loaded, rows and columns may need to be reordered. Similar for the
        # orbital coefficients.
        permutation = constructor_args.get('permutation')
        if permutation is not None:
            operators = constructor_args.get('operators')
            if operators is not None:
                for op in operators.itervalues():
                    op.apply_basis_permutation(permutation)
            wfn = constructor_args.get('wfn')
            if wfn is not None:
                wfn.expansion.apply_basis_permutation(permutation)
            del constructor_args['permutation']

        return cls(**constructor_args)

    def get_overlap(self):
        overlap = self._operators.get('olp')
        if overlap is None:
            overlap = self.lf.create_one_body(self.basis.nbasis)
            self.basis.compute_overlap(overlap)
            self._operators['olp'] = overlap
        return overlap

    def get_kinetic(self):
        kinetic = self._operators.get('kin')
        if kinetic is None:
            kinetic = self.lf.create_one_body(self.basis.nbasis)
            self.basis.compute_kinetic(kinetic)
            self._operators['kin'] = kinetic
        return kinetic

    def get_nuclear_attraction(self):

        nuclear_attraction = self._operators.get('na')
        if nuclear_attraction is None:
            nuclear_attraction = self.lf.create_one_body(self.basis.nbasis)
            # TODO: ghost atoms and extra charges
            self.basis.compute_nuclear_attraction(self.numbers.astype(float), self.coordinates, nuclear_attraction)
            self._operators['na'] = nuclear_attraction
        return nuclear_attraction

    def get_electron_repulsion(self):
        electron_repulsion = self._operators.get('er')
        if electron_repulsion is None:
            electron_repulsion = self.lf.create_two_body(self.basis.nbasis)
            self.basis.compute_electron_repulsion(electron_repulsion)
            self._operators['er'] = electron_repulsion
        return electron_repulsion
