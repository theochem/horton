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


import os
import h5py as h5, numpy as np

from horton.context import context
from horton.exceptions import ElectronCountError
from horton.grid.atgrid import AtomicGrid
from horton.grid.int1d import SimpsonIntegrator1D
from horton.grid.cext import RTransform, CubicSpline, dot_multi
from horton.guess import guess_hamiltonian_core
from horton.hamiltonian import Hamiltonian
from horton.log import log
from horton.scf import converge_scf
from horton.system import System


__all__ = ['ProAtomDB']


def _compute_average_rhos(sys, atgrid, do_population=False):
    '''Compute the spherical average of an atomic density

       **Arguments:**

       sys
            The one-atom system.

       atgrid
            The atomic grid used for computation of the spherical average.

       **Optional arguments:**

    '''
    nsphere = len(atgrid.subgrids)
    average_rhos = np.zeros(nsphere)
    for i in xrange(nsphere):
        llgrid = atgrid.subgrids[i]
        rhos = sys.compute_grid_density(llgrid.points)
        # TODO: also compute the derivatives of the average
        #       with respect to r for better splines
        average_rhos[i] = llgrid.integrate(rhos)/llgrid.weights.sum()
    if do_population:
        # TODO: ugly
        population = 4*np.pi*dot_multi(
            atgrid.rtransform.get_volume_elements(),
            atgrid.rtransform.get_radii()**2,
            average_rhos,
            atgrid.int1d.get_weights(nsphere),
        )
        return average_rhos, population
    else:
        return average_rhos


def _compute_radii(avrho, rtransform, populations):
    '''Compute approximate radii at which the atom contains the given populations

       **Arguments:**

       avrho
            A radial grid with the (spherically averaged) density of the atom.

       rtransform
            The radial transform that defines the grid.

       populations
            A list of populations for which the corresponding radii have to be
            computed.

       The return value is a list of radii corresponding to the given list of
       populations.
    '''
    # compute the integral of the density (popint)
    radii = rtransform.get_radii()
    tmp = (4*np.pi) * radii**2 * avrho * rtransform.get_volume_elements()
    popint = []
    int1d = SimpsonIntegrator1D()
    for i in xrange(len(tmp)):
        if i >= int1d.npoint_min:
            popint.append(np.dot(tmp[:i], int1d.get_weights(i)))
        else:
            popint.append(tmp[:i].sum())
    popint = np.array(popint)

    # find the radii
    indexes = popint.searchsorted(populations)
    result = []
    for i in xrange(len(populations)):
        index = indexes[i]
        if index == len(popint):
            result.append(radii[-1])
        else:
            # linear interpolation
            x = (populations[i] - popint[index])/(popint[index+1] - popint[index])
            result.append(x*radii[index+1]+(1-x)*radii[index])
    return result


class ProAtomDB(object):
    def __init__(self, rtransform, records):
        '''
           **Arguments:**

           rtransform
                An instance of one of the RTransform subclasses.

           records
                A dictionary with (number, population) keys and arrays of
                spherical averages of atomic densities.
        '''
        self._rtransform = rtransform
        self._records = records
        self._log_init()

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
        rtransform = RTransform.from_string(f.attrs['rtransform'])
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
        # TODO: add option to filter out atoms with a negative ionization potential
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
                        guess_hamiltonian_core(sys)
                    except ElectronCountError:
                        continue
                    # Compute the ground state
                    ham = Hamiltonian(sys, hamiltonian_terms, atgrid)
                    converged = converge_scf(ham)
                    if converged:
                        # Good one, store it
                        energy = ham.compute_energy()
                        average_rhos = _compute_average_rhos(sys, atgrid)
                        cases.append((energy, average_rhos))
                    elif log.do_warning:
                        log.warn('Convergence faillure when creating proatomdb from scratch. Z=%i Q=%i, M=%i' % (number, charge, mult))
                if len(cases) > 0:
                    # take the lowest in energy from all tested multiplicites
                    cases.sort()
                    records[(number, nel)] = cases[0][1]

        # Create the object.
        return cls(atgrid.rtransform, records)

    @classmethod
    def from_checkpoints(cls, fns_chk, atgrid):
        '''
           Construct a ProAtomDB from a series of Horton checkpoint files.

           **Arguments:**

           fns_chk
                A list of Horton checkpoint files.

           atgrid
                An AtomicGrid object used for computing spherical average.
        '''
        records = {}
        # Compute all caes
        for fn_chk in fns_chk:
            # Load system
            sys = System.from_file(fn_chk, chk=None) # avoid rewriting chk file
            assert sys.natom == 1
            # Compute/Get properties
            energy = sys.props.get('energy', 0.0)
            average_rhos, population = _compute_average_rhos(sys, atgrid, do_population=True)
            ipop = int(np.round(population))
            assert abs(population - ipop) < 1e-1
            key = (sys.numbers[0], ipop)
            l = records.setdefault(key, [])
            l.append((energy, average_rhos))

        # Filter out lowest energies
        for key in records.keys():
            l = records[key]
            l.sort()
            records[key] = l[0][1]

        # Create the object.
        return cls(atgrid.rtransform, records)

    @classmethod
    def from_refatoms(cls, atgrid, numbers=None, qmax=+4):
        '''
           Construct a ProAtomDB from reference atoms included in Horton

           **Arguments:**

           fns_chk
                A list of Horton checkpoint files.

           atgrid
                An AtomicGrid object used for computing spherical average.

           **Optional Arguments:**

           numbers
                A list of atom numbers to limit the selection of atoms in the
                database. When not given, all reference atoms are loaded.

           qmax
                The charge of the most positive ion allowed in the database
        '''
        # Load all the systems
        fns_chk = []
        for fn in context.glob('refatoms/*.h5'):
            name = os.path.basename(fn)
            number = int(name[:3])
            if not (numbers is None or number in numbers):
                continue
            pop = int(name[8:10])
            charge = number - pop
            if charge > qmax:
                continue
            fns_chk.append(fn)

        # Hand them over to another constructor
        return cls.from_checkpoints(fns_chk, atgrid)

    def _log_init(self):
        if log.do_medium:
            log('Initialized: %s' % self)
            log.deflist([
                ('Records', self._records.keys()),
                ('Radial Transform', self._rtransform.to_string()),
            ])
            log.blank()

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

    def get_hirshfeld_proatom_fn(self, number):
        return CubicSpline(self._records[(number,number)], rtf=self._rtransform)

    def get_hirshfeld_i_proatom_fn(self, number, pop):
        # In case of luck:
        ipop = int(pop)
        if pop == ipop:
            return CubicSpline(self._records[(number, ipop)], rtf=self._rtransform)
        else:
            del ipop

        # General case with interpolation between two integer pro-atoms.
        cpop = int(np.ceil(pop))
        fpop = int(np.floor(pop))
        if fpop == 0:
            try:
                yc = self._records[(number,cpop)]
            except KeyError:
                raise RuntimeError('No suitable proatoms found for interpolation.')
            y = pop*yc
        else:
            try:
                yc = self._records[(number,cpop)]
                yf = self._records[(number,fpop)]
            except KeyError:
                raise RuntimeError('No suitable proatoms found for interpolation.')
            y = yf*(cpop-pop) + yc*(pop-fpop)
        return CubicSpline(y, rtf=self._rtransform)

    def get_pop_range(self, number):
        pops = [p for n, p in self._records if n == number]
        return min(pops), max(pops)+1

    def compute_radii(self, number, populations, pop=None):
        '''Compute approximate radii at which the atom contains the given populations

           **Arguments:**

           number
                The element for which the radii must be computed

           populations
                A list of populations for which the corresponding radii have to
                be computed.

           **Optional argument**

           pop
                The population of the proatom. When not given, this equals
                number.

           The return value is a list of radii corresponding to the given list of
           populations.
        '''
        if pop is None:
            pop = number
        ipop = int(np.round(pop))
        avrho = self._records[(number, ipop)]
        return _compute_radii(avrho, self._rtransform, populations)
