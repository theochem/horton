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
   of a :class:`System` instance, e.g. basis sets, pseudo-potentials, ghost atoms,
   etc.
'''


import numpy as np
import h5py as h5

from horton.io import load_system_args
from horton.matrix import DenseLinalgFactory
from horton.wfn import ClosedShellWFN, OpenShellWFN


__all__ = ['System']


class System(object):
    def __init__(self, coordinates, numbers, obasis=None, wfn=None, lf=None, operators=None, props=None, chk=None):
        """
           **Arguments:**

           coordinates
                A (N, 3) float numpy array with Cartesian coordinates of the
                atoms.

           numbers
                A (N,) int numpy vector with the atomic numbers.

           obasis
                A string or an instance of either the basis set or basis set
                description classes, e.g. 'STO-3G', GOBasisDesc('STO-3G'), ...
                for the orbitals.

           **Optional arguments:**

           wfn
                A wavefunction object.

           lf
                A LinalgFactory instance. When not given, a DenseLinalgFactory
                is used by default.

           operators
                A dictionary with one- and two-body operators.

           props
                A dictionary with computed properties.

           chk
                A filename for the checkpoint file or an open h5.File object.
                If the file does not exist yet, it will be created. If the file
                already exists, it must be an HDF5 file that is structured
                such that it adheres to the format that Horton creates itself.
                If chk is an open h5.File object, it will not be closed when the
                System instance is deleted.
        """

        # A) Assign all attributes
        self._coordinates = np.array(coordinates, dtype=float, copy=False)
        self._numbers = np.array(numbers, dtype=int, copy=False)
        # some checks
        if len(self._coordinates.shape) != 2 or self._coordinates.shape[1] != 3:
            raise TypeError('coordinates argument must be a 2D array with three columns')
        if len(self._numbers.shape) != 1:
            raise TypeError('numbers must a vector of integers.')
        if self._numbers.shape[0] != self._coordinates.shape[0]:
            raise TypeError('numbers and coordinates must have compatible array shapes.')
        #
        from horton.gbasis import GOBasisDesc, GOBasis
        if isinstance(obasis, str):
            obasis_desc = GOBasisDesc(obasis)
            self._obasis = obasis_desc.apply_to(self)
        elif isinstance(obasis, GOBasisDesc):
            self._obasis = obasis.apply_to(self)
        elif isinstance(obasis, GOBasis) or obasis is None:
            self._obasis = obasis
        else:
            raise TypeError('Can not interpret %s as an orbital basis.' % obasis)
        #
        self._wfn = wfn
        #
        if lf is None:
            self._lf = DenseLinalgFactory()
        else:
            self._lf = lf
        if self._obasis is not None:
            self._lf.set_default_nbasis(self._obasis.nbasis)
        #
        if operators is None:
            self._operators = {}
        else:
            self._operators = operators
        #
        if props is None:
            self._props = {}
        else:
            self._props = props

        # Some consistency checks
        if len(self.numbers.shape) != 1 or \
           len(self.coordinates.shape) != 2 or \
           self.coordinates.shape[1] != 3 or \
           len(self.numbers) != len(self.coordinates):
            raise TypeError('Inconsistent or incorrect shapes of numbers and coordinates.')
        if self.obasis is not None:
            if self._wfn is not None and self._obasis.nbasis != self._wfn.nbasis:
                raise TypeError('The nbasis attributes of obasis and wfn are inconsistent.')
            for key, op in self._operators.iteritems():
                if op.nbasis != self._obasis.nbasis:
                    raise TypeError('The nbasis attributes of the operator %s and obasis are inconsistent.')

        # The checkpoint file
        if isinstance(chk, basestring):
            # Suppose a filename is given. Create or open an HDF5 file.
            self._chk = h5.File(chk)
            self._close_chk = True
        elif isinstance(chk, h5.File) or chk is None:
            self._chk = chk
            self._close_chk = False
        else:
            raise TypeError('The chk argument, when given, must be a filename or an open h5.File object.')
        self.update_chk()

    def __del__(self):
        # Close the HD5 checkpoint file. This must be done carefully to avoid
        # spurious error messages when an unrelated exception occurs.
        if hasattr(self, '_chk') and self.chk is not None and self._close_chk:
            self.chk.close()

    def _get_natom(self):
        '''The number of atoms'''
        return len(self.numbers)

    natom = property(_get_natom)

    def _get_coordinates(self):
        '''The positions of the nuclei'''
        return self._coordinates.view()

    coordinates = property(_get_coordinates)

    def _get_numbers(self):
        '''An array with the atomic numbers'''
        return self._numbers.view()

    numbers = property(_get_numbers)

    def _get_obasis(self):
        '''The orbital basis'''
        return self._obasis

    obasis = property(_get_obasis)

    def _get_wfn(self):
        '''The wavefunction'''
        return self._wfn

    wfn = property(_get_wfn)

    def _get_lf(self):
        '''The LinalgFactory for this system'''
        return self._lf

    lf = property(_get_lf)

    def _get_operators(self):
        '''A dictionary with operators'''
        return self._operators

    operators = property(_get_operators)

    def _get_props(self):
        '''A dictionary with computed properties.'''
        return self._props

    props = property(_get_props)

    def _get_chk(self):
        '''A ``h5py.File`` instance used as checkpoint file or ``None``'''
        return self._chk

    chk = property(_get_chk)

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

           The filenames may also contain checkpoint files and open h5.File
           objects of checkpoint files. The last such checkpoint file will
           automatically be used as a checkpoint file for this class. If you
           want to override this behavior, provide the ``chk`` keyword argument
           (may be None).
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
                wfn.apply_basis_permutation(permutation)
            del constructor_args['permutation']

        return cls(**constructor_args)

    def update_chk(self, field_name=None):
        """Write (a part of) the system to the checkpoint file.

           **Optional Argument:**

           field
                A field string that specifies which part must be written to the
                checkpoint file. When not given, all possible fields are
                written. The latter is only useful in specific cases, e.g. upon
                initialization of the system. The available field names are
                specified in the attribute register dictionary in the
                module ``horton.checkpoint``.
        """
        if self._chk is not None:
            from horton.checkpoint import register
            if field_name is None:
                for field_name, field in register.iteritems():
                    field.write(self._chk, self)
            else:
                field = register[field_name]
                field.write(self._chk, self)

    def init_wfn(self, charge=0, mult=None, restricted=None):
        '''Initialize a wavefunction object.

           **Optional Arguments:**

           charge
                The total charge of the system. Defaults to zero.

           mult
                The spin multiplicity. Defaults to lowest possible.

           restricted
                Set to True or False to enforce a restricted or unrestricted
                wavefunction. Note that restricted open shell is not yet
                supported. When not set, restricted is used when mult==1 and
                unrestricted otherwise.
        '''
        if self._wfn is not None:
            raise RuntimeError('A wavefunction is already present.')
        if self._obasis is None:
            raise RuntimeError('A wavefunction can only be initialized when a basis is specified.')
        if charge is None:
            charge = 0
        nel = self.numbers.sum() - charge
        if mult is None:
            mult = nel%2+1
        elif ((nel%2 == 0) ^ (mult%2 != 0)):
            raise ValueError('Not compatible: number of electrons = %i and spin multiplicity = %i' % (nel, mult))
        if mult < 1:
            raise ValueError('mult must be strictly positive.')
        if restricted is True and mult != 1:
            raise ValueError('Restricted==True only works when mult==1. Restricted open shell is not supported yet.')
        if restricted is None:
            restricted = mult==1
        if restricted:
            self._wfn = ClosedShellWFN(nel/2, self.lf, self.obasis.nbasis)
        else:
            nalpha = (nel + (mult-1))/2
            nbeta = (nel - (mult-1))/2
            self._wfn = OpenShellWFN(nalpha, nbeta, self.lf, self.obasis.nbasis)

    def get_overlap(self):
        overlap = self._operators.get('olp')
        if overlap is None:
            overlap = self.lf.create_one_body()
            self.obasis.compute_overlap(overlap)
            self._operators['olp'] = overlap
            self.update_chk('operators.olp')
        return overlap

    def get_kinetic(self):
        kinetic = self._operators.get('kin')
        if kinetic is None:
            kinetic = self.lf.create_one_body()
            self.obasis.compute_kinetic(kinetic)
            self._operators['kin'] = kinetic
            self.update_chk('operators.kin')
        return kinetic

    def get_nuclear_attraction(self):
        nuclear_attraction = self._operators.get('na')
        if nuclear_attraction is None:
            nuclear_attraction = self.lf.create_one_body()
            # TODO: ghost atoms and extra charges
            self.obasis.compute_nuclear_attraction(self.numbers.astype(float), self.coordinates, nuclear_attraction)
            self._operators['na'] = nuclear_attraction
            self.update_chk('operators.na')
        return nuclear_attraction

    def get_electron_repulsion(self):
        electron_repulsion = self._operators.get('er')
        if electron_repulsion is None:
            electron_repulsion = self.lf.create_two_body()
            self.obasis.compute_electron_repulsion(electron_repulsion)
            self._operators['er'] = electron_repulsion
            # ER integrals are not checkpointed by default because they are too heavy.
            # Can be done manually by user if needed: ``system.update_chk('operators.er')``
            #self.update_chk('operators.er')
        return electron_repulsion

    def compute_density_grid(self, points, use_dm=False, rhos=None, select='full'):
        '''Compute the electron density on a grid

           **Arguments:**

           points
                A Numpy array with grid points, shape (npoint,3)

           **Optional arguments:**

           use_dm
                Use C routine based on the density matrix instead of the orbitals.

           rhos
                An output array, shape (npoint,). The results are added to this
                array.

           select
                'alpha', 'beta', 'full' or 'spin'. ('full' is the default.)

           **Returns:**

           rhos
                The array with the densities. This is the same as the output
                argument, in case it was provided.
        '''
        if rhos is None:
            rhos = np.zeros(len(points), float)
        elif rhos.shape != (points.shape[0],):
            raise TypeError('The shape of the output array is wrong')
        if use_dm:
            dm = self.lf.create_one_body()
            self.wfn.compute_density_matrix(dm, select=select)
            self.obasis.compute_density_grid_dm(dm, points, rhos)
        else:
            for expansion, nocc, scale in self.wfn.iter_expansions(select):
                self.obasis.compute_density_grid_orb(expansion, nocc, scale, points, rhos)
        return rhos

    def compute_nucnuc(self):
        '''Compute interaction energy of the nuclei'''
        # TODO: move this to low-level code one day.
        result = 0.0
        for i in xrange(self.natom):
            for j in xrange(i):
                distance = np.linalg.norm(self.coordinates[i]-self.coordinates[j])
                result += self.numbers[i]*self.numbers[j]/distance
        self._props['energy_nn'] = result
        return result

    def compute_grid_one_body(self, points, weights, pots, one_body):
        '''See documentation self.obasis.compute_grid_one_body'''
        self.obasis.compute_grid_one_body(points, weights, pots, one_body)
