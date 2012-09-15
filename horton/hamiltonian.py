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


import numpy as np

from horton.cache import Cache


__all__ = [
    'Hamiltonian', 'HamiltonianTerm', 'KineticEnergy', 'ExternalPotential',
    'Hartree', 'HartreeFock', 'DiracExchange'
]


class Hamiltonian(object):
    def __init__(self, system, terms, grid=None):
        '''
           **Arguments:**

           system
                The System object for which the energy must be computed.

           terms
                The terms in the Hamiltonian. Kinetic energy and external
                potential (nuclei) are added automatically.

           **Optional arguments:**

           grid
                The integration grid, in case some terms need one.
        '''
        # check arguments:
        if len(terms) == 0:
            raise ValueError('At least one term must be present in the Hamiltonian.')
        for term in terms:
            if term.require_grid and grid is None:
                raise TypeError('The term %s requires a grid, but not grid is given.' % term)

        # Assign attributes
        self.system = system
        self.terms = list(terms)
        self.grid = grid

        # Add standard terms if missing
        #  1) Kinetic energy
        if sum(isinstance(term, KineticEnergy) for term in terms) == 0:
            self.terms.append(KineticEnergy())
        #  2) Hartree (or HatreeFock, which is a subclass of Hartree)
        if sum(isinstance(term, Hartree) for term in terms) == 0:
            self.terms.append(Hartree())
        #  3) External Potential
        if sum(isinstance(term, ExternalPotential) for term in terms) == 0:
            self.terms.append(ExternalPotential())

        # Create a cache for shared intermediate results.
        self.cache = Cache()

        # Pre-compute stuff
        for term in self.terms:
            term.prepare_system(self.system, self.cache, self.grid)

        # Compute overlap matrix
        self.overlap = system.get_overlap()

    def invalidate_derived(self):
        '''Mark the properties derived from the wfn as outdated.

           This method does not recompute anything, but just marks operators
           and density matrices as outdated. They are recomputed as they are
           needed.
        '''
        self.cache.invalidate()
        for dm in self.system.dms.itervalues():
            dm.invalidate()
        for term in self.terms:
            term.invalidate_derived()

    def compute_energy(self):
        '''Compute energy.

           **Returns:**

           The total energy, including nuclear-nuclear repulsion.
        '''
        total = 0.0
        self.system.update_dms()
        dm_alpha = self.system.dms.get('alpha')
        dm_beta = self.system.dms.get('beta')
        dm_full = self.system.dms.get('full')
        for term in self.terms:
            total += term.compute_energy(dm_alpha, dm_beta, dm_full)
        total += self.system.compute_nucnuc()
        # Store result in chk file
        self.system._props['energy'] = total
        self.system.update_chk('props')
        return total

    def compute_fock(self, fock_alpha, fock_beta):
        '''Compute alpha (and beta) Fock matrix(es).

           **Arguments:**

           fock_alpha
                A One-Body operator output argument for the alpha fock matrix.

           fock_alpha
                A One-Body operator output argument for the beta fock matrix.

           The density matrices in ``self.system.dms`` are used for the
           computations, so make sure they are up to date, i.e. call
           ``system.update_dms()`` prior to ``compute_energy``.

           In the case of a closed-shell computation, the argument fock_beta is
           ``None``.
        '''
        self.system.update_dms()
        dm_alpha = self.system.dms.get('alpha')
        dm_beta = self.system.dms.get('beta')
        dm_full = self.system.dms.get('full')
        for term in self.terms:
            term.add_fock_matrix(dm_alpha, dm_beta, dm_full, fock_alpha, fock_beta)


class HamiltonianTerm(object):
    require_grid = False
    def prepare_system(self, system, cache, grid):
        self.system = system
        self.cache = cache
        self.grid = grid

    def invalidate_derived(self):
        pass

    # Generic update routines that may be useful to various base classes
    def update_rho(self, select):
        rho, new = self.cache.load('rho_%s' % select, alloc=self.grid.size)
        if new:
            self.system.compute_density_grid(self.grid.points, rhos=rho, select=select)
        return rho


    def compute_energy(self, dm_alpha, dm_beta, dm_full):
        raise NotImplementedError

    def add_fock_matrix(self, dm_alpha, dm_beta, dm_full, fock_alpha, fock_beta):
        raise NotImplementedError


class FixedTerm(HamiltonianTerm):
    '''Base class for all terms that are linear in the density matrix

       This is (technically) a special class because the Fock operator does not
       have to be recomputed when the density matrix changes.
    '''
    def get_operator(self, system):
        # subclasses should return the operator and a suffix for the energy_* property.
        raise NotImplementedError

    def prepare_system(self, system, cache, grid):
        HamiltonianTerm.prepare_system(self, system, cache, grid)
        self.operator, self.suffix = self.get_operator(system)

    def compute_energy(self, dm_alpha, dm_beta, dm_full):
        if dm_beta is None:
            result = 2*self.operator.expectation_value(dm_alpha)
        else:
            result = self.operator.expectation_value(dm_full)
        self.system._props['energy_%s' % self.suffix] = result
        return result

    def add_fock_matrix(self, dm_alpha, dm_beta, dm_full, fock_alpha, fock_beta):
        for fock in fock_alpha, fock_beta:
            if fock is not None:
                fock.iadd(self.operator, 1)


class KineticEnergy(FixedTerm):
    def get_operator(self, system):
        return system.get_kinetic(), 'kin'


class ExternalPotential(FixedTerm):
    def get_operator(self, system):
        tmp = system.get_nuclear_attraction()
        tmp.iscale(-1)
        return tmp, 'ne'


class Hartree(HamiltonianTerm):
    def prepare_system(self, system, cache, grid):
        HamiltonianTerm.prepare_system(self, system, cache, grid)
        self.electron_repulsion = system.get_electron_repulsion()
        self.coulomb = system.lf.create_one_body(system.obasis.nbasis)

    def invalidate_derived(self):
        self.coulomb.invalidate()

    def _update_coulomb(self, dm_alpha, dm_beta, dm_full):
        '''Recompute the Coulomb operator if it has become invalid'''
        if not self.coulomb.valid:
            if dm_beta is None:
                self.electron_repulsion.apply_direct(dm_alpha, self.coulomb)
                self.coulomb.iscale(2)
            else:
                self.electron_repulsion.apply_direct(dm_full, self.coulomb)

    def compute_energy(self, dm_alpha, dm_beta, dm_full):
        self._update_coulomb(dm_alpha, dm_beta, dm_full)
        if dm_beta is None:
            result = self.coulomb.expectation_value(dm_alpha)
        else:
            result = 0.5*self.coulomb.expectation_value(dm_full)
        self.system._props['energy_hartree'] = result
        return result

    def add_fock_matrix(self, dm_alpha, dm_beta, dm_full, fock_alpha, fock_beta):
        self._update_coulomb(dm_alpha, dm_beta, dm_full)
        if dm_beta is None:
            # closed shell
            fock_alpha.iadd(self.coulomb, 1)
        else:
            # open shell
            fock_alpha.iadd(self.coulomb, 1)
            fock_beta.iadd(self.coulomb, 1)


class HartreeFock(Hartree):
    def __init__(self, fraction_exchange=1.0):
        self.fraction_exchange = fraction_exchange

    def prepare_system(self, system, cache, grid):
        Hartree.prepare_system(self, system, cache, grid)
        self.exchange_alpha = system.lf.create_one_body(system.obasis.nbasis)
        if not system.wfn.closed_shell:
            self.exchange_beta = system.lf.create_one_body(system.obasis.nbasis)
        else:
            self.exchange_beta = None

    def invalidate_derived(self):
        Hartree.invalidate_derived(self)
        self.exchange_alpha.invalidate()
        if self.exchange_beta is not None:
            self.exchange_beta.invalidate()

    def _update_exchange(self, dm_alpha, dm_beta, dm_full):
        '''Recompute the Exchange operator(s) if invalid'''
        if not self.exchange_alpha.valid:
            self.electron_repulsion.apply_exchange(dm_alpha, self.exchange_alpha)
        if self.exchange_beta is not None and not self.exchange_beta.valid:
            self.electron_repulsion.apply_exchange(dm_beta, self.exchange_beta)

    def compute_energy(self, dm_alpha, dm_beta, dm_full):
        energy_hartree = Hartree.compute_energy(self, dm_alpha, dm_beta, dm_full)
        self._update_exchange(dm_alpha, dm_beta, dm_full)
        if dm_beta is None:
            energy_fock = -self.exchange_alpha.expectation_value(dm_alpha)
        else:
            energy_fock = -0.5*self.exchange_alpha.expectation_value(dm_alpha) \
                          -0.5*self.exchange_beta.expectation_value(dm_beta)
        self.system._props['energy_exchange_fock'] = energy_fock
        return energy_hartree + self.fraction_exchange*energy_fock

    def add_fock_matrix(self, dm_alpha, dm_beta, dm_full, fock_alpha, fock_beta):
        Hartree.add_fock_matrix(self, dm_alpha, dm_beta, dm_full, fock_alpha, fock_beta)
        self._update_exchange(dm_alpha, dm_beta, dm_full)
        fock_alpha.iadd(self.exchange_alpha, -self.fraction_exchange)
        if fock_beta is not None:
            fock_beta.iadd(self.exchange_beta, -self.fraction_exchange)


class DiracExchange(HamiltonianTerm):
    '''A (for now naive) implementation of the Dirac Exchange Functional'''

    require_grid = True
    def __init__(self, coeff=None):
        '''
           **Arguments:**



           **Optional arguments:**

           coeff
                The coefficient Cx in front of the Dirac exchange energy.
                It defaults to the uniform electron gas value, i.e.
                Cx = 3/4 (3/pi)^(1/3).
        '''
        if coeff is None:
            self.coeff = 3.0/4.0*(3.0/np.pi)**(1.0/3.0)
        else:
            self.coeff = coeff
        self.derived_coeff = -self.coeff*(4.0/3.0)*2**(1.0/3.0)

    def prepare_system(self, system, cache, grid):
        HamiltonianTerm.prepare_system(self, system, cache, grid)
        self.exchange = {}
        self.exchange['alpha'] = system.lf.create_one_body(system.obasis.nbasis)
        if not self.system.wfn.closed_shell:
            self.exchange['beta'] = system.lf.create_one_body(system.obasis.nbasis)

    def invalidate_derived(self):
        for exchange in self.exchange.itervalues():
            exchange.invalidate()

    def _update_exchange(self, dm_alpha, dm_beta, dm_full):
        '''Recompute the Exchange operator(s) if invalid'''
        def helper(select):
            # update grid stuff
            rho = self.update_rho(select)
            pot, new = self.cache.load('pot_exchange_dirac_%s' % select, alloc=self.grid.size)
            if new:
                pot[:] = self.derived_coeff*rho**(1.0/3.0)

            # update operator stuff
            # TODO: move operators to cache.
            exchange = self.exchange[select]
            if not exchange.valid:
                exchange.reset()
                self.system.compute_grid_one_body(self.grid.points, self.grid.weights, pot, exchange)

        helper('alpha')
        if not self.system.wfn.closed_shell:
            helper('beta')

    def compute_energy(self, dm_alpha, dm_beta, dm_full):
        self._update_exchange(dm_alpha, dm_beta, dm_full)

        def helper(select):
            pot = self.cache.load('pot_exchange_dirac_%s' % select)
            rho = self.cache.load('rho_%s' % select)
            return self.grid.integrate(pot, rho)

        energy = helper('alpha')
        if not self.system.wfn.closed_shell:
            energy += helper('beta')
        else:
            energy *= 2
        energy *= 3.0/4.0
        self.system._props['energy_exchange_dirac'] = energy
        return energy

    def add_fock_matrix(self, dm_alpha, dm_beta, dm_full, fock_alpha, fock_beta):
        self._update_exchange(dm_alpha, dm_beta, dm_full)
        fock_alpha.iadd(self.exchange['alpha'])
        if fock_beta is not None:
            fock_beta.iadd(self.exchange['beta'])
