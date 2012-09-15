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
'''Tools for playing with pro-atomic databases'''


import h5py as h5, numpy as np

from horton.exceptions import ElectronCountError
from horton.grid.atgrid import AtomicGrid
from horton.grid.rtransform import BaseRTransform
from horton.guess import guess_hamiltonian_core
from horton.hamiltonian import Hamiltonian
from horton.scf import converge_scf
from horton.system import System


__all__ = ['ProAtomDB']


class ProAtomDB(object):
    def __init__(self, rtransform, records):
        '''
           **Arguments:**

           rtransform
                An instance of one of the BaseRTransform subclasses.

           records
                A dictionary with (number, population) keys and arrays of
                spherical averages of atomic densities.
        '''
        self._rtransform = rtransform
        self._records = records

    @classmethod
    def from_file(cls, filename):
        '''Construct an dabase from an HDF5 file

           **Arguments:**

           filename
                A string with the filename of the hdf5 file, or a h5.File or
                h5.Group object.
        '''
        # parse the argument
        if isinstance(filename, basestring):
            f = h5.File(filename, 'r')
            do_close = True
        elif isinstance(filename, h5.Group) or isinstance(filename, h5.File):
            f = filename
            do_close = False
        # Read
        rtransform = BaseRTransform.from_string(f.attrs['rtransform'], None)
        records = {}
        for name, dataset in f.iteritems():
            number, population = name.split(':')
            number = int(number)
            population = int(population)
            records[(number, population)] = np.array(dataset)
        return ProAtomDB(rtransform, records)
        # close
        if do_close:
            f.close()

    @classmethod
    def from_scratch(cls, hamiltonian_terms, obasis, atgrid, numbers, qmin=-2, qmax=3, mult_range=6):
        '''
           **Arguments:**

           hamiltonian_terms
                A list of HamiltonianTerm instances.

           obasis
                The orbital basis to use for the atomic computations.

           atgrid
                An AtomicGrid object used for integration.

           numbers
                A list of atomic numbers for the database.

           **Optional arguments:**

           qmin
                The most negative charge to consider.

           qmax
                The most positive charge to consider.

           mult_range
                The maximum multiplicity when scanning for the multiplicity
                with the lowest energy. Hund's rule is not used because it
                may not be valid for ions and because we are lazy!
        '''
        # TODO: allow for different grids for each element.
        # TODO long term: make it possible to control the scf code and to update
        #                 an existing db with additional comps.
        # check the essentials of the arguments:
        if not isinstance(atgrid, AtomicGrid) or atgrid.subgrids is None:
            raise TypeError('The subgrids of the atomic grid are needed. Use the argument keep_subgrids=1 when creating the atomic grid.')

        # Compute all records, considering all reasonable multiplicites
        records = {}
        for number in numbers:
            for charge in xrange(qmin, qmax+1):
                # select the range of multiplicities to consider
                nel = number-charge
                if nel%2 == 0:
                    mults = xrange(1, (mult_range/2+1)*2, 2)
                else:
                    mults = xrange(2, (mult_range/2)*2+1, 2)
                # compute each spin state
                cases = []
                for mult in mults:
                    # setup the system ...
                    sys = System(np.array([atgrid.center]), np.array([number]), obasis=obasis)
                    # ... and the initial wavefn
                    try:
                        sys.init_wfn(charge, mult)
                    except ElectronCountError:
                        continue
                    guess_hamiltonian_core(sys)
                    # Compute the ground state
                    ham = Hamiltonian(sys, hamiltonian_terms, atgrid)
                    converged = converge_scf(ham)
                    if converged:
                        # Good one, store it
                        energy = ham.compute_energy()
                        nsphere = len(atgrid.subgrids)
                        average_rhos = np.zeros(nsphere)
                        for i in xrange(nsphere):
                            llgrid = atgrid.subgrids[i]
                            rhos = sys.compute_density_grid(llgrid.points)
                            # TODO: also compute the derivatives of the average
                            #       with respect to r for better splines
                            average_rhos[i] = llgrid.integrate(rhos)/llgrid.weights.sum()
                        cases.append((energy, average_rhos))
                    # TODO: when a screenlogger is implemented, it must show
                    # a warning when one of the cases does not converge.
                if len(cases) > 0:
                    # take the lowest in energy from all tested multiplicites
                    cases.sort()
                    records[(number, nel)] = cases[0][1]

        # Create the object.
        return cls(atgrid.rtransform, records)

    def to_file(self, filename):
        '''Write the database to an HDF5 file


           **Arguments:**

           filename
                A string with the filename of the hdf5 file, or a h5.File or
                h5.Group object.
        '''
        # parse the argument

        if isinstance(filename, basestring):
            f = h5.File(filename, 'w')
            do_close = True
        elif isinstance(filename, h5.Group) or isinstance(filename, h5.File):
            f = filename
            do_close = False
        # Write
        f.attrs['rtransform'] = self._rtransform.to_string()
        for key, average_rho in self._records.iteritems():
            f.create_dataset('%i:%i' % key, data=average_rho)
        # close
        if do_close:
            f.close()
