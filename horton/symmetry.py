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
'''Geometric symmetries

   The symmetry tools in HORTON are just meant to provide optional additional
   information on top of a System instance.
'''


import numpy as np

from horton.cext import Cell
from horton.exceptions import SymmetryError
from horton.periodic import periodic


__all__ = ['Symmetry']


class Symmetry(object):
    '''An geometrical symmetry descriptor'''
    def __init__(self, name, generators, fracs, numbers, cell, labels=None):
        '''
           **Arguments:**

           name
                Whatever name you want to give to this symmetry. This is
                converted to a string.

           generators
                A list of (3,4) matrices where the first three columns contain
                the linear transformation matrix and the last column is a
                translation vector. These transformations act on the fractional
                (or reduced) coordinates.

           fracs
                The fractional coordinates of a primitive cell/unit.

           numbers
                The corresponding element numbers of the primitive cell/unit

           cell
                A Cell object. Even for isolated systems a cell object must
                be provided with nvec=0.

           **Optional arguments:**

           labels
                A list of unique labels for each atom in the primitive unit.
                These are generated from the numbers array when not given.
        '''
        for generator in generators:
            if not (isinstance(generator, np.ndarray) and generator.shape == (3, 4)):
                raise TypeError('the list of generators may only contain 3x4 arrays.')
        if not (isinstance(fracs, np.ndarray) and fracs.ndim == 2 and fracs.shape[1] == 3):
            raise TypeError('fracs must be a numpy array with three columns.')
        if not (isinstance(numbers, np.ndarray) and numbers.ndim == 1):
            raise TypeError('numbers must be a numpy one-dimensional array.')
        if numbers.shape[0] != fracs.shape[0]:
            raise TypeError('numbers and fracs are not consistent in size.')
        if not isinstance(cell, Cell):
            raise TypeError('The cell object must be of the Cell class.')
        if labels is None:
            labels = []
            for i, n in enumerate(numbers):
                labels.append('%s%i' % (periodic[n].symbol, i))
        else:
            if len(labels) != numbers.shape[0]:
                raise TypeError('The length of the labels list is not consistent with the number of atoms in the primitive unit.')
            if len(set(labels)) < len(labels):
                raise ValueError('Not all labels are unique.')

        self._name = str(name)
        self._generators = np.array(generators)
        self._fracs = fracs
        self._numbers = numbers
        self._cell = cell
        self._labels = np.array(labels)

    @classmethod
    def from_hdf5(cls, grp):
        '''Construct a Symmetry object from a HDF5 group'''
        return cls(
            grp['name'][()],
            grp['generators'][:],
            grp['fracs'][:],
            grp['numbers'][:],
            Cell.from_hdf5(grp['cell']),
            grp['labels'][:],
        )

    def to_hdf5(self, grp):
        '''Write a symmetry object to an HDF5 group'''
        grp.attrs['class'] = self.__class__.__name__
        grp['name'] = self.name
        grp['generators'] = self.generators
        grp['fracs'] = self.fracs
        grp['numbers'] = self.numbers
        tmp = grp.create_group('cell')
        self.cell.to_hdf5(tmp)
        grp['labels'] = self.labels


    def _get_name(self):
        return self._name

    name = property(_get_name)

    def _get_generators(self):
        return self._generators

    generators = property(_get_generators)

    def _get_natom(self):
        '''The number of atoms in the primitive unit'''
        return self.numbers.shape[0]

    natom = property(_get_natom)

    def _get_fracs(self):
        return self._fracs

    fracs = property(_get_fracs)

    def _get_numbers(self):
        return self._numbers

    numbers = property(_get_numbers)

    def _get_cell(self):
        return self._cell

    cell = property(_get_cell)

    def _get_labels(self):
        return self._labels

    labels = property(_get_labels)

    def _iter_images(self, cell=None):
        '''Loops over all periodic images of all atoms in the primitive cell

           This loop may include atoms whose symmetric images coincide.
        '''
        if cell is None:
            cell = self.cell
        for i in xrange(self.natom):
            for j in xrange(len(self.generators)):
                g = self.generators[j]
                frac = np.dot(g[:,:3], self.fracs[i]) + g[:,3]
                cart = cell.to_cart(frac)
                yield i, j, g, frac, cart


    def generate(self, threshold=0.001):
        '''Apply the generators to the primitive unit to obtain a full molecule.

           **Optional arguments:**

           threshold
                When, after transformation with the generators, two (or more)
                symmetry atoms overlap within this threshold, they will be
                merged into one. The distance is measured in Cartesian
                coordinates.

           **Returns:**

           coordinates
                Cartesian coordinates for all atoms.

           numbers
                Atomic numbers for all atoms.

           links
                An array of indexes to connect each atom back with an atom in
                the primitive cell (first column) and a generator (second
                column).
        '''
        coordinates = []
        numbers = []
        links = []
        for i, j, g, frac, cart in self._iter_images():
            # test if it is already present
            duplicate = False
            for ocart in coordinates:
                delta = ocart - cart
                self.cell.mic(delta)
                if np.linalg.norm(delta) < threshold:
                    duplicate = True
                    break
            if duplicate:
                continue

            coordinates.append(cart)
            numbers.append(self.numbers[i])
            links.append([i, j])

        return np.array(coordinates), np.array(numbers), np.array(links)

    def identify(self, coordinates, cell, threshold=0.1):
        '''Connect atoms in the primitive unit with atoms in the full molecule

           **Arguments:**

           coordinates
                An (N, 3) array of atomic coordinates that adhere (with some
                minor deviation) to this symmetry.

           cell
                A Cell instance describing the periodic boundary conditions

           **Optional arguments:**

           threshold
                The maximum allowed distance between the ideal atom position
                and the actual atom position

           **Returns:**

           links
                An array of indexes to connect each atom back with an atom in
                the primitive cell (first column) and a generator (second
                column).

           If an atom in the full molecule can not be linked with an atom in
           the primitive unit, a SymmetryError is raised. If the full molecule
           contains less atoms than the perfect crystal (e.g. a vacancy), this
           method will not complain.
        '''
        if len(coordinates.shape) != 2 or coordinates.shape[1] != 3:
            raise TypeError('The argument coordinates must be an array with three columns.')
        natom = coordinates.shape[0]
        links = []
        for k in xrange(natom):
            match = False
            for i, j, g, frac, cart in self._iter_images(cell):
                delta = coordinates[k] - cart
                cell.mic(delta)
                #print i, j, np.linalg.norm(delta)
                if np.linalg.norm(delta) < threshold:
                    match = True
                    break
            if not match:
                raise SymmetryError('Could not find atom in primitive unit corresponding to atom %i' % k)

            links.append([i, j])
        return np.array(links)
