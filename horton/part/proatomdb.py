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
'''Pro-atom databases'''


import os
import h5py as h5, numpy as np

from horton.context import context
from horton.grid.atgrid import AtomicGrid, AtomicGridSpec
from horton.grid.cext import RTransform, CubicSpline
from horton.grid.radial import RadialGrid
from horton.io.lockedh5 import LockedH5File
from horton.log import log, timer
from horton.io.iodata import IOData


__all__ = ['ProAtomRecord', 'ProAtomDB']


class ProAtomRecord(object):
    '''A single proatomic density record'''

    @classmethod
    def from_iodata(cls, iodata, agspec='fine'):
        '''Construct a proatom record from a data dictionary and an atomic grid

           **Arguments:**

           iodata
                IOData instance

           **Optional arguments:**

           agspec
                A specifications of the atomic grid. This can either be an
                instance of the AtomicGridSpec object, or the first argument
                of its constructor.
        '''
        if iodata.natom != 1:
            raise ValueError('More than one atom found for pro-atom record.')
        number = iodata.numbers[0]
        pseudo_number = iodata.pseudo_numbers[0]
        center = iodata.coordinates[0]
        energy = getattr(iodata, 'energy', 0.0)
        dm_full = iodata.get_dm_full()
        return cls.from_dm(center, number, pseudo_number, iodata.obasis, dm_full, energy, agspec)

    @classmethod
    def from_dm(cls, center, number, pseudo_number, obasis, dm_full, energy, agspec='fine'):
        '''Construct a proatom record from a single-atom density matrix

           **Arguments:**

           center
                The position of the nuclues (needed for spherical average)

           number
                The atomic number

           pseudo_number
                The effective core charges

           obasis
                The orbital basis

           dm_full
                The spin-summed density matrix object

           energy
                The electronic energy of the atom

           **Optional arguments:**

           agspec
                A specifications of the atomic grid. This can either be an
                instance of the AtomicGridSpec object, or the first argument
                of its constructor.
        '''
        if len(center) != 3:
            raise TypeError('Center should be a vector with three elements.')
        if not isinstance(agspec, AtomicGridSpec):
            agspec = AtomicGridSpec(agspec)

        # Compute the density and the gradient on a grid
        atgrid = AtomicGrid(number, pseudo_number, center, agspec)
        rho_all = obasis.compute_grid_density_dm(dm_full, atgrid.points)
        grad_all = obasis.compute_grid_gradient_dm(dm_full, atgrid.points)
        rho, deriv = atgrid.get_spherical_average(rho_all, grads=[grad_all])

        # Derive the number of electrions and the charge
        overlap = dm_full.new()
        obasis.compute_overlap(overlap)
        nel = overlap.contract_two('ab,ba', dm_full)
        assert abs(nel - int(np.round(nel))) < 1e-4 # only integer nel are supported
        nel = int(np.round(nel))
        charge = pseudo_number - nel

        # Create object
        return cls(number, charge, energy, atgrid.rgrid, rho, deriv, pseudo_number)

    def __init__(self, number, charge, energy, rgrid, rho, deriv=None, pseudo_number=None, ipot_energy=None):
        '''
           **Arguments:**

           number
                The element number of the proatom.

           charge
                The net charge of the proatom. (integer)

           energy
                The total energy of the proatom.

           rgrid
                The radial grid on which the density (and optional radial
                density derivatives) are tabulated.

           rho
                The electron density on the grid.

           **Optional arguments:**

           deriv
                The radial derivative of the electron density.

           pseudo_number
                The effective core charge (defaults to number).

           ipot_energy
                The ionization potential.
        '''
        self._number = number
        self._charge = charge
        self._energy = energy
        self._rho = rho
        self._deriv = deriv
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

    def _get_ipot_energy(self):
        '''The ionization potential'''
        return self._ipot_energy

    ipot_energy = property(_get_ipot_energy)

    def _get_rho(self):
        '''The density on a radial grid'''
        return self._rho

    rho = property(_get_rho)

    def _get_deriv(self):
        '''The radial derivative of the density on a radial grid'''
        return self._deriv

    deriv = property(_get_deriv)

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
        tmp = (4*np.pi) * radii**2 * self.rho * self.rgrid.rtransform.get_deriv()
        popint = []
        for i in xrange(len(tmp)):
            if i >= self.rgrid.int1d.npoint_min:
                popint.append(np.dot(tmp[:i], self.rgrid.int1d.get_weights(i)))
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
        if self._deriv is not None:
            self._deriv = self._deriv[:npoint]
        self._rgrid = self._rgrid.chop(npoint)

    def __eq__(self, other):
        return (self.number == other.number and
                self.charge == other.charge and
                self.energy == other.energy and
                (self.rho == other.rho).all() and
                ((self.deriv is None and other.deriv is None) or (self.deriv is not None and other.deriv is not None and self.deriv == other.deriv).all()) and
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
        # Search for duplicates (same number and charge) and only retain the
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
    def from_files(cls, fns, agspec='fine'):
        '''
           Construct a ProAtomDB from a series of HORTON checkpoint files.

           **Arguments:**

           fns_chk
                A list of atomic output files.

           **Optional arguments:**

           agspec
                A specifications of the atomic grid. This can either be an
                instance of the AtomicGridSpec object, or the first argument
                of its constructor.
        '''
        if not isinstance(agspec, AtomicGridSpec):
            agspec = AtomicGridSpec(agspec)
        records = []
        for fn in fns:
            # Load atomic data
            with timer.section('Load proatom'):
                mol = IOData.from_file(fn)
            with timer.section('Proatom grid'):
                records.append(ProAtomRecord.from_iodata(mol, agspec))
        return cls(records)

    @classmethod
    def from_refatoms(cls, numbers=None, max_cation=3, max_anion=2, agspec='fine'):
        '''
           Construct a ProAtomDB from reference atoms included in HORTON

           **Arguments:**

           fns_chk
                A list of HORTON checkpoint files.

           **Optional Arguments:**

           numbers
                A list of atom numbers to limit the selection of atoms in the
                database. When not given, all reference atoms are loaded.

           max_cation
                The charge of the most positive cation to include

           max_anion
                Minus the charge of the most negativ anion to include. (This
                is typically a positive argument.)

           agspec
                A specifications of the atomic grid. This can either be an
                instance of the AtomicGridSpec object, or the first argument
                of its constructor.
        '''
        # Search for all the relevant .h5 files of built-in reference atoms
        fns_chk = []
        for fn in context.glob('refatoms/*.h5'):
            name = os.path.basename(fn)
            number = int(name[:3])
            if not (numbers is None or number in numbers):
                continue
            pop = int(name[8:10])
            charge = number - pop
            if charge > max_cation or charge < -max_anion:
                continue
            fns_chk.append(fn)

        # Hand them over to another constructor
        return cls.from_files(fns_chk, agspec)

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
        if isinstance(filename, h5.Group):
            records = load_proatom_records_h5_group(filename)
        elif isinstance(filename, basestring):
            if filename.endswith('.h5'):
                records = load_proatom_records_h5_file(filename)
            elif filename.endswith('.atdens'):
                records = load_proatom_records_atdens(filename)
            else:
                raise ValueError('Proatomdb file type not supported')
        else:
            raise NotImplementedError

        return ProAtomDB(records)

    def to_file(self, filename):
        '''Write the database to an HDF5 file


           **Arguments:**

           filename
                A string with the filename of the hdf5 file, or a h5.File or
                h5.Group object.
        '''
        # parse the argument
        if isinstance(filename, basestring):
            f = LockedH5File(filename, 'w')
            do_close = True
        elif isinstance(filename, h5.Group):
            f = filename
            do_close = False
        try:
            # Write
            for record in self._records:
                name = 'Z=%i_Q=%+i' % (record.number, record.charge)
                if name in f:
                    del f[name]
                grp = f.create_group(name)
                grp.attrs['number'] = record.number
                grp.attrs['charge'] = record.charge
                grp.attrs['energy'] = record.energy
                grp.attrs['rtransform'] = record.rgrid.rtransform.to_string()
                grp['rho'] = record.rho
                if record.deriv is not None:
                    grp['deriv'] = record.deriv
                grp.attrs['pseudo_number'] = record.pseudo_number
                if record.ipot_energy is not None:
                    grp.attrs['ipot_energy'] = record.ipot_energy
        finally:
            # close
            if do_close:
                f.close()

    def get_rho(self, number, parameters=0, combine='linear', do_deriv=False):
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

           do_deriv
                When set to True, the derivative of rho is also returned. In
                case the derivative is not available, the second return value is
                None.
        '''
        if isinstance(parameters, int):
            charge = parameters
            record = self.get_record(number, charge)
            if do_deriv:
                return record.rho, record.deriv
            else:
                return record.rho
        elif isinstance(parameters, dict):
            rho = 0.0
            if do_deriv:
                deriv = 0.0
            else:
                deriv = None
            if combine == 'linear':
                for charge, coeff in parameters.iteritems():
                    if coeff != 0.0:
                        record = self.get_record(number, charge)
                        rho += coeff*record.rho
                        if do_deriv and record.deriv is not None and deriv is not None:
                            deriv += coeff*record.deriv
                        else:
                            deriv = None
            elif combine == 'geometric':
                for charge, coeff in parameters.iteritems():
                    if coeff != 0.0:
                        record = self.get_record(number, charge)
                        rho += coeff*np.log(record.rho)
                        if do_deriv and record.deriv is not None and deriv is not None:
                            deriv += coeff*record.deriv/record.rho
                        else:
                            deriv = None
                rho = np.exp(rho)
                if do_deriv and deriv is not None:
                    deriv = rho*deriv
            else:
                raise ValueError('Combine argument "%s" not supported.' % combine)
            if not isinstance(rho, np.ndarray):
                rho = self.get_rgrid(number).zeros()
                if do_deriv:
                    deriv = self.get_rgrid(number).zeros()
            if do_deriv:
                return rho, deriv
            else:
                return rho
        else:
            raise TypeError('Could not interpret parameters argument')

    def get_spline(self, number, parameters=0, combine='linear'):
        '''Construct a proatom spline.

           **Arguments:** See ``get_rho`` method.
        '''
        rho, deriv = self.get_rho(number, parameters, combine, do_deriv=True)
        return CubicSpline(rho, deriv, self.get_rgrid(number).rtransform)

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


def load_proatom_records_h5_group(f):
    '''Load proatom records from the given HDF5 group'''
    records = []
    for grp in f.itervalues():
        assert isinstance(grp, h5.Group)
        if 'deriv' in grp:
            deriv = grp['deriv'][:]
        else:
            deriv = None
        records.append(ProAtomRecord(
            number=grp.attrs['number'],
            charge=grp.attrs['charge'],
            energy=grp.attrs['energy'],
            rgrid=RadialGrid(RTransform.from_string(grp.attrs['rtransform'])),
            rho=grp['rho'][:],
            deriv=deriv,
            pseudo_number=grp.attrs.get('pseudo_number'),
            ipot_energy=grp.attrs.get('ipot_energy'),
        ))
    return records


def load_proatom_records_h5_file(filename):
    '''Load proatom records from the given HDF5 file'''
    with LockedH5File(filename) as f:
        return load_proatom_records_h5_group(f)


def load_proatom_records_atdens(filename):
    '''Load proatom records from the given atdens file file'''
    def read_numbers(f, npoint):
        numbers = []
        while len(numbers) < npoint:
            words = f.next().replace('D', 'E').split()
            for word in words:
                numbers.append(float(word))
        return np.array(numbers)

    from horton.grid.cext import ExpRTransform
    records = []
    with open(filename) as f:
        # load the radii
        words = f.next().split()
        assert words[0] == 'RADII'
        npoint = int(words[1])
        radii = read_numbers(f, npoint)
        # Construct a radial grid
        r1 = radii[1]
        r2 = radii[-1]
        rtf = ExpRTransform(r1, r2, npoint-1)
        rgrid = RadialGrid(rtf)
        # load the proatoms
        while True:
            try:
                words = f.next().split()
                number = int(words[0])
                population = int(words[1])
                charge = number - population
                rho = read_numbers(f, npoint)[1:]
                record = ProAtomRecord(number, charge, 0.0, rgrid, rho)
                records.append(record)
            except StopIteration:
                break
    return records
