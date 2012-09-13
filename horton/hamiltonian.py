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

    def compute_energy(self, dm_alpha, dm_beta, dm_full):
        '''Compute energy.

           **Arguments:**

           dm_alpha
                The density matrix of the spin-up electrons

           dm_beta
                The density matrix of the spin-down electrons

           dm_full
                The density matrix of all the electrons

           **Returns:**

           The total energy, including nuclear-nuclear repulsion.

           In the case of a closed-shell computation, the arguments dm_beta,
           dm_full are ``None``.
        '''
        # TODO: store all sorts of energies in system object and checkpoint file
        total = 0.0
        for term in self.terms:
            total += term.compute_energy(dm_alpha, dm_beta, dm_full)
        total += self.system.compute_nucnuc()
        return total

    def compute_fock(self, dm_alpha, dm_beta, dm_full, fock_alpha, fock_beta):
        '''Compute alpha (and beta) Fock matrix(es).

           **Arguments:**

           dm_alpha
                The density matrix of the spin-up electrons

           dm_beta
                The density matrix of the spin-down electrons

           dm_full
                The density matrix of all the electrons

           fock_alpha
                A One-Body operator output argument for the alpha fock matrix.

           fock_alpha
                A One-Body operator output argument for the beta fock matrix.

           In the case of a closed-shell computation, the arguments dm_beta,
           dm_full and fock_beta are ``None``.
        '''
        for term in self.terms:
            term.add_fock_matrix(dm_alpha, dm_beta, dm_full, fock_alpha, fock_beta)


class HamiltonianTerm(object):
    def prepare_system(self, system):
        pass

    def compute_energy(self, dm_alpha, dm_beta, dm_full):
        raise NotImplementedError

    def add_fock_matrix(self, dm_alpha, dm_beta, dm_full, fock_alpha, fock_beta):
        raise NotImplementedError


class KineticEnergy(HamiltonianTerm):
    def prepare_system(self, system):
        self.kinetic = system.get_kinetic()

    def compute_energy(self, dm_alpha, dm_beta, dm_full):
        if dm_beta is None:
            return 2*self.kinetic.expectation_value(dm_alpha)
        else:
            return self.kinetic.expectation_value(dm_full)

    def add_fock_matrix(self, dm_alpha, dm_beta, dm_full, fock_alpha, fock_beta):
        for fock in fock_alpha, fock_beta:
            if fock is not None:
                fock.iadd(self.kinetic, 1)


class Hartree(HamiltonianTerm):
    def prepare_system(self, system):
        self.electron_repulsion = system.get_electron_repulsion()
        self.coulomb = system.lf.create_one_body(system.obasis.nbasis)

    def compute_energy(self, dm_alpha, dm_beta, dm_full):
        if dm_beta is None:
            self.electron_repulsion.apply_direct(dm_alpha, self.coulomb)
            self.coulomb.iscale(2)
            return self.coulomb.expectation_value(dm_alpha)
        else:
            self.electron_repulsion.apply_direct(dm_full, self.coulomb)
            return 0.5*self.coulomb.expectation_value(dm_full)

    def add_fock_matrix(self, dm_alpha, dm_beta, dm_full, fock_alpha, fock_beta):
        if dm_beta is None:
            # closed shell
            self.electron_repulsion.apply_direct(dm_alpha, self.coulomb)
            self.coulomb.iscale(2)
            fock_alpha.iadd(self.coulomb, 1)
        else:
            # open shell
            self.electron_repulsion.apply_direct(dm_full, self.coulomb)
            fock_alpha.iadd(self.coulomb, 1)
            fock_beta.iadd(self.coulomb, 1)


class HartreeFock(Hartree):
    def __init__(self, fraction_exchange=1.0):
        self.fraction_exchange = fraction_exchange

    def prepare_system(self, system):
        Hartree.prepare_system(self, system)
        self.exchange = system.lf.create_one_body(system.obasis.nbasis)

    def compute_energy(self, dm_alpha, dm_beta, dm_full):
        result = Hartree.compute_energy(self, dm_alpha, dm_beta, dm_full)
        if dm_beta is None:
            self.electron_repulsion.apply_exchange(dm_alpha, self.exchange)
            return result - self.fraction_exchange*self.exchange.expectation_value(dm_alpha)
        else:
            self.electron_repulsion.apply_exchange(dm_alpha, self.exchange)
            result -= 0.5*self.fraction_exchange*self.exchange.expectation_value(dm_alpha)
            self.electron_repulsion.apply_exchange(dm_beta, self.exchange)
            result -= 0.5*self.fraction_exchange*self.exchange.expectation_value(dm_beta)
            return result

    def add_fock_matrix(self, dm_alpha, dm_beta, dm_full, fock_alpha, fock_beta):
        Hartree.add_fock_matrix(self, dm_alpha, dm_beta, dm_full, fock_alpha, fock_beta)
        for dm, fock in (dm_alpha, fock_alpha), (dm_beta, fock_beta):
            if dm is not None:
                self.electron_repulsion.apply_exchange(dm, self.exchange)
                fock.iadd(self.exchange, -self.fraction_exchange)


class ExternalPotential(HamiltonianTerm):
    def prepare_system(self, system):
        self.nuclear_attraction = system.get_nuclear_attraction()

    def compute_energy(self, dm_alpha, dm_beta, dm_full):
        if dm_beta is None:
            return -2*self.nuclear_attraction.expectation_value(dm_alpha)
        else:
            return -self.nuclear_attraction.expectation_value(dm_full)

    def add_fock_matrix(self, dm_alpha, dm_beta, dm_full, fock_alpha, fock_beta):
        for fock in fock_alpha, fock_beta:
            if fock is not None:
                fock.iadd(self.nuclear_attraction, -1)
