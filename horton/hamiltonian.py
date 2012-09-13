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


__all__ = [
    'Hamiltonian', 'HamiltonianTerm', 'KineticEnergy', 'Hartree', 'HartreeFock',
    'ExternalPotential',
]


class Hamiltonian(object):
    def __init__(self, system, terms):
        self.system = system
        self.terms = list(terms)
        if len(terms) == 0:
            raise ValueError('At least one term must be present in the Hamiltonian.')
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
        # Pre-compute stuff
        for term in self.terms:
            term.prepare_system(self.system)
        # Compute overlap matrix
        self.overlap = system.get_overlap()

    def invalidate_derived(self):
        '''Mark the properties derived from the wfn as outdated.

           This method does not recompute anything, but just marks operators
           and density matrices as outdated. They are recomputed as they are
           needed.
        '''
        for dm in self.system.dms.itervalues():
            dm.invalidate()
        for term in self.terms:
            term.invalidate_derived()

    def compute_energy(self):
        '''Compute energy.

           **Returns:**

           The total energy, including nuclear-nuclear repulsion.
        '''
        # TODO: store all sorts of energies in system object and checkpoint file
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
    def prepare_system(self, system):
        self.system = system

    def invalidate_derived(self):
        pass

    def compute_energy(self, dm_alpha, dm_beta, dm_full):
        raise NotImplementedError

    def add_fock_matrix(self, dm_alpha, dm_beta, dm_full, fock_alpha, fock_beta):
        raise NotImplementedError


class KineticEnergy(HamiltonianTerm):
    def prepare_system(self, system):
        HamiltonianTerm.prepare_system(self, system)
        self.kinetic = system.get_kinetic()

    def compute_energy(self, dm_alpha, dm_beta, dm_full):
        if dm_beta is None:
            result = 2*self.kinetic.expectation_value(dm_alpha)
        else:
            result = self.kinetic.expectation_value(dm_full)
        self.system._props['energy_kin'] = result
        return result


    def add_fock_matrix(self, dm_alpha, dm_beta, dm_full, fock_alpha, fock_beta):
        for fock in fock_alpha, fock_beta:
            if fock is not None:
                fock.iadd(self.kinetic, 1)


class Hartree(HamiltonianTerm):
    def prepare_system(self, system):
        HamiltonianTerm.prepare_system(self, system)
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

    def prepare_system(self, system):
        Hartree.prepare_system(self, system)
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


# TODO: derive from common base class with kinetic energy
class ExternalPotential(HamiltonianTerm):
    def prepare_system(self, system):
        HamiltonianTerm.prepare_system(self, system)
        self.nuclear_attraction = system.get_nuclear_attraction()

    def compute_energy(self, dm_alpha, dm_beta, dm_full):
        if dm_beta is None:
            result = -2*self.nuclear_attraction.expectation_value(dm_alpha)
        else:
            result = -self.nuclear_attraction.expectation_value(dm_full)
        self.system._props['energy_ne'] = result
        return result

    def add_fock_matrix(self, dm_alpha, dm_beta, dm_full, fock_alpha, fock_beta):
        for fock in fock_alpha, fock_beta:
            if fock is not None:
                fock.iadd(self.nuclear_attraction, -1)
