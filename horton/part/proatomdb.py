# -*- coding: utf-8 -*-
# Horton is a Density Functional Theory program.
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
'''Tools for playing with pro-atomic databases'''


import os
import h5py as h5, numpy as np

from horton.context import context
from horton.exceptions import ElectronCountError
from horton.grid.atgrid import AtomicGrid
from horton.grid.int1d import SimpsonIntegrator1D
from horton.grid.cext import RTransform, CubicSpline, dot_multi
from horton.grid.radial import RadialGrid
from horton.guess import guess_hamiltonian_core
from horton.hamiltonian import Hamiltonian
from horton.log import log, timer
from horton.scf import converge_scf
from horton.system import System


__all__ = ['ProAtomRecord', 'ProAtomDB']


class ProAtomRecord(object):
    @classmethod
    def from_system(cls, system, atgrid):
        '''Construct a proatom record from a system and an atomic grid

           **Arguments:**

           sys
                The one-atom system.

           atgrid
                The atomic grid used for computation of the spherical average.
        '''
        if system.natom != 1:
            raise ValueError('Can only construct a proatom record from a one-atom system.')

        # Compute the spherically averaged density
        rho_full = system.compute_grid_density(atgrid.points)
        rho = atgrid.get_spherical_average(rho_full)

        # Derive the number of electrions
        try:
            nel = system.wfn.nel
        except NotImplementedError:
            nel = int(np.round(atgrid.rgrid.integrate(rho)))

        # Get all the other information from the atom
        energy = system.props['energy']
        try:
            homo_energy = system.wfn.homo_energy
        except AttributeError:
            homo_energy = None
        number = system.numbers[0]
        pseudo_number = system.pseudo_numbers[0]
        charge = pseudo_number - nel

        # Create object
        return cls(number, charge, energy, homo_energy, atgrid.rgrid, rho, pseudo_number)

    def __init__(self, number, charge, energy, homo_energy, rgrid, rho, pseudo_number=None, ipot_energy=None):
        self._number = number
        self._charge = charge
        self._energy = energy
        self._homo_energy = homo_energy
        self._rho = rho
        self._rgrid = rgrid
        if pseudo_number is None:
            self._pseudo_number = number
        else:
            self._pseudo_number = pseudo_number
        if self.pseudo_population == 1:
            self._ipot_energy = -self._energy
        else:
            self._ipot_energy = ipot_energy
        self._safe = True

    def _get_number(self):
        '''The element number'''
        return self._number

    number = property(_get_number)

    def _get_charge(self):
        '''The charge'''
        return self._charge

    charge = property(_get_charge)

    def _get_energy(self):
        '''The total electronic energy'''
        return self._energy

    energy = property(_get_energy)

    def _get_homo_energy(self):
        '''The HOMO energy'''
        return self._homo_energy

    homo_energy = property(_get_homo_energy)

    def _get_ipot_energy(self):
        '''The ionization potential'''
        return self._ipot_energy

    ipot_energy = property(_get_ipot_energy)

    def _get_rho(self):
        '''The density on a radial grid'''
        return self._rho

    rho = property(_get_rho)

    def _get_rgrid(self):
        '''The radial grid'''
        return self._rgrid

    rgrid = property(_get_rgrid)

    def _get_pseudo_number(self):
        '''The pseudo element number (effective core charge)'''
        return self._pseudo_number

    pseudo_number = property(_get_pseudo_number)

    def _get_population(self):
        '''The total number of electrons'''
        return self._number - self._charge

    population = property(_get_population)

    def _get_pseudo_population(self):
        '''The total effective number of electrons'''
        return self._pseudo_number - self._charge

    pseudo_population = property(_get_pseudo_population)

    def _get_safe(self):
        '''When safe is True, this pro atom is safe to use, i.e. not know to be basis-set bound'''
        return self._safe

    safe = property(_get_safe)

    def update_safe(self, other):
        '''Updates the safe attribute based on a comparison with other records

           **Arguments:**

           other
                Another instance of ProAtomRecord
        '''
        if other.number == self._number:
            if other.population < self.population and \
               other.energy < self.energy:
                self._safe = False
            if other.population == self.population - 1:
                self._ipot_energy = other.energy - self.energy

    def compute_radii(self, populations):
        '''Compute approximate radii and grid points at which the atom contains the given populations

           **Arguments:**

           populations
                A list of populations for which the corresponding radii have to
                be computed.
        '''
        # compute the running integral of the density (popint)
        radii = self.rgrid.radii
        tmp = (4*np.pi) * radii**2 * self.rho * self.rgrid.rtransform.get_volume_elements()
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
                x = (populations[i] - popint[index])/(popint[index-1] - popint[index])
                result.append(x*radii[index-1]+(1-x)*radii[index])
        return indexes, result

    def get_moment(self, order):
        '''Return the integral of rho*r**order'''
        return self.rgrid.integrate(self.rho, self.rgrid.radii**order)

    def chop(self, npoint):
        '''Reduce the proatom to the given number of radial grid points.'''
        self._rho = self._rho[:npoint]
        self._rgrid = self._rgrid.chop(npoint)

    def __eq__(self, other):
        return (self.number == other.number and
                self.charge == other.charge and
                self.energy == other.energy and
                self.homo_energy == other.homo_energy and
                (self.rho == other.rho).all() and
                self.rgrid == other.rgrid and
                self.pseudo_number == other.pseudo_number and
                self.ipot_energy == other.ipot_energy)

    def __ne__(self, other):
        return not self.__eq__(other)


class ProAtomDB(object):
    def __init__(self, records):
        '''
           **Arguments:**

           records
                A list of ProAtomRecord instances. If two or more records have
                the same number and charge, only the lowest in energy is
                retained.

           Based on the records present it is determined which records are
           safe to use, i.e. apparently not bound by the basis set.
        '''
        # Search for diplicates (same number and charge) and only retain the
        # lowest in energy for each combination.
        _map = {}
        for r in records:
            l = _map.setdefault((r.number, r.charge), [])
            l.append(r)
        for key, l in _map.iteritems():
            l.sort(key=(lambda r: r.energy))
        records = [l[0] for l in _map.itervalues()]

        # Store attribtues
        self._records = records
        self._map = dict(((r.number, r.charge), r) for r in records)

        # check that all records of a given element have the same rgrid
        self._rgrid_map = {}
        for r in records:
            rgrid = self._rgrid_map.get(r.number)
            if rgrid is None:
                # store
                self._rgrid_map[r.number] = r.rgrid
            else:
                # compare
                if rgrid != r.rgrid:
                    raise ValueError('All proatoms of a given element must have the same radial grid')

        # Update the safe flags based on the energies of other pro_atoms
        for number in self.get_numbers():
            for charge0 in self.get_charges(number):
                r0 = self.get_record(number, charge0)
                for charge1 in self.get_charges(number):
                    r1 = self.get_record(number, charge1)
                    r0.update_safe(r1)

        # Screen info
        self._log_init()

    def _log_init(self):
        if log.do_medium:
            log('Initialized: %s' % self)
            log.deflist([
                ('Numbers', self._rgrid_map.keys()),
                ('Records', self._map.keys()),
            ])
            log.blank()

    def get_record(self, number, charge):
        return self._map[(number, charge)]

    def get_numbers(self):
        '''Return the element numbers present in the database'''
        result = self._rgrid_map.keys()
        result.sort()
        return result

    def get_charges(self, number, safe=False):
        result = [r.charge for r in self._records if r.number == number and (r.safe or (not safe))]
        result.sort(reverse=True)
        return result

    def get_rgrid(self, number):
        return self._rgrid_map[number]

    def _get_size(self):
        return len(self._records)

    size = property(_get_size)

    @classmethod
    def from_files(cls, fns, atgrid):
        '''
           Construct a ProAtomDB from a series of Horton checkpoint files.

           **Arguments:**

           fns_chk
                A list of system files.

           atgrid
                An AtomicGrid object used for computing spherical average.
        '''
        # TODO: add option to use standard grids
        records = []
        for fn in fns:
            # Load system
            with timer.section('Load proatom'):
                sys = System.from_file(fn, chk=None) # avoid rewriting in case of chk file
            with timer.section('Proatom grid'):
                records.append(ProAtomRecord.from_system(sys, atgrid))
        return cls(records)

    @classmethod
    def from_refatoms(cls, atgrid, numbers=None, max_kation=3, max_anion=2):
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

           max_kation
                The charge of the most positive kation to include

           max_anion
                Minus the charge of the most negativ anion to include. (This
                is typically a positive argument.)
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
            if charge > max_kation or charge < -max_anion:
                continue
            fns_chk.append(fn)

        # Hand them over to another constructor
        return cls.from_files(fns_chk, atgrid)

    @classmethod
    def from_file(cls, filename):
        '''Construct an dabase from an HDF5 file

           **Arguments:**

           filename
                A string with the filename of the hdf5 file, or a h5.File or
                h5.Group object.

           Note that the records are loaded and given as argument to the
           constructor, which may weed out duplicates.
        '''
        # parse the argument
        if isinstance(filename, basestring):
            f = h5.File(filename, 'r')
            do_close = True
        elif isinstance(filename, h5.Group):
            f = filename
            do_close = False
        # Read
        records = []
        for grp in f.itervalues():
            assert isinstance(grp, h5.Group)
            records.append(ProAtomRecord(
                number=grp.attrs['number'],
                charge=grp.attrs['charge'],
                energy=grp.attrs['energy'],
                homo_energy=grp.attrs.get('homo_energy'),
                rgrid=RadialGrid(RTransform.from_string(grp.attrs['rtransform'])),
                rho=grp['rho'][:],
                pseudo_number=grp.attrs.get('pseudo_number'),
                ipot_energy=grp.attrs.get('ipot_energy'),
            ))
        result = ProAtomDB(records)
        # close
        if do_close:
            f.close()
        return result

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
        elif isinstance(filename, h5.Group):
            f = filename
            do_close = False
        # Write
        for record in self._records:
            name = 'Z=%i_Q=%+i' % (record.number, record.charge)
            if name in f:
                del f[name]
            grp = f.create_group(name)
            grp.attrs['number'] = record.number
            grp.attrs['charge'] = record.charge
            grp.attrs['energy'] = record.energy
            if record.homo_energy is not None:
                grp.attrs['homo_energy'] = record.homo_energy
            grp.attrs['rtransform'] = record.rgrid.rtransform.to_string()
            grp['rho'] = record.rho
            grp.attrs['pseudo_number'] = record.pseudo_number
            if record.ipot_energy is not None:
                grp.attrs['ipot_energy'] = record.ipot_energy
        # close
        if do_close:
            f.close()

    def get_rho(self, number, parameters=0, combine='linear'):
        '''Construct a proatom density on a grid.

           **Arguments:**

           number
                The element

           **Optional arguments:**

           parameters
                This argument may take two forms:

                * Integer: the charge of the pro-atom

                * Dictionary: a linear or geometric combination of different
                  charged pro-atoms. The keys are the charges and the values
                  are the coefficients.

           combine
                In case parameters is an array, this determines the type of
                combination: 'linear' or 'geometric'.
        '''
        if isinstance(parameters, int):
            charge = parameters
            record = self.get_record(number, charge)
            return record.rho
        elif isinstance(parameters, dict):
            rho = 0.0
            if combine == 'linear':
                for charge, coeff in parameters.iteritems():
                    if coeff != 0.0:
                        rho += coeff*self.get_record(number,charge).rho
            elif combine == 'geometric':
                for charge, coeff in parameters.iteritems():
                    if coeff != 0.0:
                        rho += coeff*np.log(self.get_record(number,charge).rho)
                rho = np.exp(rho)
            else:
                raise ValueError('Combine argument "%s" not supported.' % combine)
            if not isinstance(rho, np.ndarray):
                rho = self.get_rgrid(number).zeros()
            return rho
        else:
            raise TypeError('Could not interpret parameters argument')


    def get_spline(self, number, parameters=0, combine='linear'):
        '''Construct a proatom spline.

           **Arguments:** See ``get_rho.. method.
        '''
        rho = self.get_rho(number, parameters, combine)
        return CubicSpline(rho, rtf=self.get_rgrid(number).rtransform)

    def compact(self, nel_lost):
        '''Make the pro-atoms more compact

           **Argument:**

           nel_lost
                This parameter controls the part of the tail that gets
                neglected. This is the (maximym) number of electrons in the part
                that gets discarded.

           Note that only 'safe' atoms are considered to determine the cutoff
           radius.
        '''
        if log.do_medium:
            log('Reducing extents of the pro-atoms')
            log('   Z     npiont           radius')
            log.hline()
        for number in self.get_numbers():
            rgrid = self.get_rgrid(number)
            npoint = 0
            for charge in self.get_charges(number, safe=True):
                r = self.get_record(number, charge)
                nel = r.pseudo_number-charge
                npoint = max(npoint, r.compute_radii([nel-nel_lost])[0][0]+1)
            for charge in self.get_charges(number):
                r = self.get_record(number, charge)
                r.chop(npoint)
            new_rgrid = self._rgrid_map[number].chop(npoint)
            self._rgrid_map[number] = new_rgrid
            if log.do_medium:
                log('%4i   %5i -> %5i    %10.3e -> %10.3e' % (
                    number, rgrid.size, new_rgrid.size, rgrid.radii[-1],
                    new_rgrid.radii[-1]
                ))
        if log.do_medium:
            log.hline()

    def normalize(self):
        if log.do_medium:
            log('Normalizing proatoms to integer populations')
            log('   Z  charge             before             after')
            log.hline()
        for number in self.get_numbers():
            rgrid = self.get_rgrid(number)
            for charge in self.get_charges(number):
                r = self.get_record(number, charge)
                nel_before = rgrid.integrate(r.rho)
                nel_integer = r.pseudo_number - charge
                r.rho[:] *= nel_integer/nel_before
                nel_after = rgrid.integrate(r.rho)
                if log.do_medium:
                    log('%4i     %+3i    %15.8e   %15.8e' % (
                        number, charge, nel_before, nel_after
                    ))
