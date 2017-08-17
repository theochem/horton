# -*- coding: utf-8 -*-
# HORTON: Helpful Open-source Research TOol for N-fermion systems.
# Copyright (C) 2011-2017 The HORTON Development Team
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
"""Pro-atom databases"""


from __future__ import print_function

import os
import h5py as h5
import numpy as np

from horton.grid import RTransform, CubicSpline


__all__ = ['ProAtomRecord', 'ProAtomDB']


class ProAtomRecord(object):
    """A single proatomic density record"""

    def __init__(self, number, charge, energy, rgrid, rho, deriv=None, pseudo_number=None, ipot_energy=None):
        """
        Parameters
        ----------
        number : int
            The atomic number of the proatom.
        charge : int
            The net charge of the proatom.
        energy : float
            The total energy of the proatom.
        rgrid : instance RadialGrid
            The radial grid on which the density (and optional radial
            density derivatives) are tabulated.
        rho : np.ndarray
            The electron density on the grid.
        deriv : np.ndarray
            The radial derivative of the electron density.
        pseudo_number : int, default=None
            The effective core charge. If None, it is set to number.
        ipot_energy : float, default=None
            The ionization potential.
        """
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

    @property
    def number(self):
        """The atomic number."""
        return self._number

    @property
    def pseudo_number(self):
        """The pseudo atomic number (effective core charge)"""
        return self._pseudo_number

    @property
    def charge(self):
        """The atomic charge."""
        return self._charge

    @property
    def energy(self):
        """The total electronic energy."""
        return self._energy

    @property
    def ipot_energy(self):
        """The ionization potential."""
        return self._ipot_energy

    @property
    def rho(self):
        """The density on a radial grid."""
        return self._rho

    @property
    def deriv(self):
        """The radial derivative of the density on a radial grid."""
        return self._deriv

    @property
    def rgrid(self):
        """The radial grid."""
        return self._rgrid

    @property
    def population(self):
        """The total number of electrons."""
        return self._number - self._charge

    @property
    def pseudo_population(self):
        """The total effective number of electrons."""
        return self._pseudo_number - self._charge

    @property
    def safe(self):
        """When safe is True, this pro atom is safe to use, i.e. not know to be basis-set bound"""
        return self._safe

    def update_safe(self, other):
        """Update the safe attribute based on a comparison with other records.

           **Arguments:**

           other
                Another instance of ProAtomRecord
        """
        if other.number == self._number:
            if other.population < self.population and \
               other.energy < self.energy:
                self._safe = False
            if other.population == self.population - 1:
                self._ipot_energy = other.energy - self.energy

    def compute_radii(self, populations):
        """Compute approximate radii and grid points at which the atom contains the given populations

           **Arguments:**

           populations
                A list of populations for which the corresponding radii have to
                be computed.
        """
        # compute the running integral of the density (popint)
        radii = self.rgrid.radii
        tmp = 4 * np.pi * radii**2 * self.rho * self.rgrid.rtransform.get_deriv()
        popint = tmp.cumsum()
        # find the radii
        indexes = popint.searchsorted(populations)
        result = []
        for i in range(len(populations)):
            index = indexes[i]
            if index == len(popint):
                result.append(radii[-1])
            else:
                # linear interpolation
                x = (populations[i] - popint[index]) / (popint[index - 1] - popint[index])
                result.append(x * radii[index - 1] + (1 - x) * radii[index])
        return indexes, result

    def get_moment(self, order):
        """Return the integral of rho*r**order"""
        return self.rgrid.integrate(self.rho, self.rgrid.radii**order)

    def chop(self, npoint):
        """Reduce the proatom to the given number of radial grid points."""
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
        """
        Parameters
        ----------
        records : sequence of ProAtomRecord instances
            A sequence of ProAtomRecord instances. If two or more records have
            the same number and charge, only the lowest in energy is
            retained.

           Based on the records present it is determined which records are
           safe to use, i.e. apparently not bound by the basis set.
        """
        # Search for duplicates (same number and charge) and only retain the
        # lowest in energy for each combination.
        _map = {}
        for r in records:
            l = _map.setdefault((r.number, r.charge), [])
            l.append(r)
        for key, l in list(_map.items()):
            l.sort(key=(lambda r: r.energy))
        records = [l[0] for l in list(_map.values())]

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
                    raise ValueError('All proatoms of a given element must have the same radial '
                                     'grid.')

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
        print('5: Initialized: %s' % self)
        print([
            ('5: Numbers', list(self._rgrid_map.keys())),
            ('5: Records', list(self._map.keys())),
        ])
        print()

    def get_record(self, number, charge):
        return self._map[(number, charge)]

    def get_numbers(self):
        """Return the element numbers present in the database"""
        result = list(self._rgrid_map.keys())
        result.sort()
        return result

    def get_charges(self, number, safe=False):
        result = [r.charge for r in self._records if r.number == number and (r.safe or (not safe))]
        result.sort(reverse=True)
        return result

    def get_rgrid(self, number):
        return self._rgrid_map[number]

    @property
    def size(self):
        """Number of proatoms in the database."""
        return len(self._records)

    def get_rho(self, number, parameters=0, combine='linear', do_deriv=False):
        """Construct a proatom density on a grid.

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
        """
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
                for charge, coeff in list(parameters.items()):
                    if coeff != 0.0:
                        record = self.get_record(number, charge)
                        rho += coeff * record.rho
                        if do_deriv and record.deriv is not None and deriv is not None:
                            deriv += coeff * record.deriv
                        else:
                            deriv = None
            elif combine == 'geometric':
                for charge, coeff in list(parameters.items()):
                    if coeff != 0.0:
                        record = self.get_record(number, charge)
                        rho += coeff * np.log(record.rho)
                        if do_deriv and record.deriv is not None and deriv is not None:
                            deriv += coeff * record.deriv / record.rho
                        else:
                            deriv = None
                rho = np.exp(rho)
                if do_deriv and deriv is not None:
                    deriv = rho * deriv
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
        """Construct a proatom spline.

           **Arguments:** See ``get_rho`` method.
        """
        rho, deriv = self.get_rho(number, parameters, combine, do_deriv=True)
        return CubicSpline(rho, deriv, self.get_rgrid(number).rtransform)

    def compact(self, nel_lost):
        """Make the pro-atoms more compact

           **Argument:**

           nel_lost
                This parameter controls the part of the tail that gets
                neglected. This is the (maximym) number of electrons in the part
                that gets discarded.

           Note that only 'safe' atoms are considered to determine the cutoff
           radius.
        """
        print('5:Reducing extents of the pro-atoms')
        print('5:   Z     npiont           radius')
        for number in self.get_numbers():
            rgrid = self.get_rgrid(number)
            npoint = 0
            for charge in self.get_charges(number, safe=True):
                r = self.get_record(number, charge)
                nel = r.pseudo_number - charge
                npoint = max(npoint, r.compute_radii([nel - nel_lost])[0][0] + 1)
            for charge in self.get_charges(number):
                r = self.get_record(number, charge)
                r.chop(npoint)
            new_rgrid = self._rgrid_map[number].chop(npoint)
            self._rgrid_map[number] = new_rgrid
            print('5:%4i   %5i -> %5i    %10.3e -> %10.3e' % (
                number, rgrid.size, new_rgrid.size, rgrid.radii[-1], new_rgrid.radii[-1]))
        print()

    def normalize(self):
        print('5:Normalizing proatoms to integer populations')
        print('5:   Z  charge             before             after')
        print()
        for number in self.get_numbers():
            rgrid = self.get_rgrid(number)
            for charge in self.get_charges(number):
                r = self.get_record(number, charge)
                nel_before = rgrid.integrate(r.rho)
                nel_integer = r.pseudo_number - charge
                r.rho[:] *= nel_integer / nel_before
                nel_after = rgrid.integrate(r.rho)
                print('5:%4i     %+3i    %15.8e   %15.8e' % (number, charge, nel_before, nel_after))
