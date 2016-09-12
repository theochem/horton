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
"""Mean-field DFT/HF Hamiltonian data structures"""


from horton.log import log
from horton.cache import Cache
from horton.utils import doc_inherit


__all__ = [
    'REffHam', 'UEffHam'
]


class EffHam(object):
    """Base class for the effective Hamiltonians.

    Class attributes
    ----------------
    ndm : int
        The number of input density matrices and output fock matrices (e.g. ndm=1 for
        restricted wfns, ndm=2 for unrestricted wfns.)
    deriv_scale : float
        In principle, the fock matrix is the derivative of the expectation value towards
        the density matrix elements. In practice, this is not always the case. Depending
        on the type of effective Hamiltonian, the fock matrices must be multiplied with a
        factor to obtain proper derivatives. This factor is stored in the class attribute
        ``deriv_scale``. It defaults to 1.0.
    """

    ndm = None
    deriv_scale = 1.0

    def __init__(self, terms, external=None):
        """Initialize an EffHam instance.

        Parameters
        ----------
        terms : list with instances of Observable
            The terms in the Hamiltonian.
        external : dict
            A dictionary with external energy contributions that do not depend on the
            wavefunction, e.g. nuclear-nuclear interactions or QM/MM mechanical embedding
            terms. Use ``nn`` as key for the nuclear-nuclear term.
        """
        # check arguments:
        if len(terms) == 0:
            raise ValueError('At least one term must be present in the Hamiltonian.')

        # Assign attributes
        self.terms = list(terms)
        self.external = {} if external is None else external

        # Create a cache for shared intermediate results. This cache should only
        # be used for derived quantities that depend on the wavefunction and
        # need to be updated at each SCF cycle.
        self.cache = Cache()

    def reset(self, *dms):
        """Remove intermediate results from cache and specify new input density matrices.

        Parameters
        ----------
        dm1, dm2, ... : TwoIndex
            The input density matrices. Their interpretation is fixed in
            derived classes.
        """
        raise NotImplementedError

    def reset_delta(self, *delta_dms):
        """Remove intermediate results for delta_dms from cache and specify new inputs.

        Parameters
        ----------
        delta_dm1, delta_dm2, ... : TwoIndex
            First-order changes to the density matrix, used to compute the dot product
            with the energy Hessian.
        """
        raise NotImplementedError

    def compute_energy(self):
        """Compute the total energy.

        The input for this method must be provided through the ``reset``
        method.

        Returns
        -------
        energy : float
            The expectation value, including the constant terms defined through the
            ``external`` argument of the constructor
        """
        total = 0.0
        for term in self.terms:
            energy = term.compute_energy(self.cache)
            self.cache['energy_%s' % term.label] = energy
            total += energy
        for key, energy in self.external.iteritems():
            self.cache['energy_%s' % key] = energy
            total += energy
        self.cache['energy'] = total
        return total

    def log(self):
        """Write an overview of the last computation on screen."""
        log('Contributions to the energy:')
        log.hline()
        log('                                              term                 Value')
        log.hline()
        for term in self.terms:
            energy = self.cache['energy_%s' % term.label]
            log('%50s  %20.12f' % (term.label, energy))
        for key, energy in self.external.iteritems():
            log('%50s  %20.12f' % (key, energy))
        log('%50s  %20.12f' % ('total', self.cache['energy']))
        log.hline()
        log.blank()

    def compute_fock(self, *focks):
        """Compute the fock matrices.

        A Fock matrix is the derivative of the energy toward the components of the
        corresponding input density matrices.

        The input for this method must be provided through the ``reset`` method. Note that
        the Fock matrix must be multiplied by the factor deriv_scale to obtain the proper
        derivative of the energy toward the density matrix, in order to maintain the
        common conventions for Fock matrices.

        Parameters
        ----------
        fock1, fock2, ... : TwoIndex
            A list of output Fock operators. Old content is discarded.
        """
        for fock in focks:
            fock.clear()
        # Loop over all terms and add contributions to the Fock matrix.
        for term in self.terms:
            term.add_fock(self.cache, *focks)

    def compute_dot_hessian(self, *outputs):
        """Compute the dot product of the energy Hessian with a delta DM.

        The Hessian in this method is the second derivative of the energy towards the
        matrix elements of the input density matrix or matrices. The ``dms`` and
        ``delta_dms`` are set via the ``reset`` and ``reset_delta`` methods, respectively.

        Parameters
        ----------
        outputs : TwoIndex
            A list of output TwoIndex objects in which the dot product of the energy
            Hessian with the delta density matrices is stored.

        Note that the result must be multiplied by the feactor deriv_scale squared in
        order to obtain the proper second order derivative. This is due to conventions
        related to the definition of the Fock matrix.
        """
        for output in outputs:
            output.clear()
        # Loop over all terms and add contributions to the output two-index
        # objects.
        for term in self.terms:
            term.add_dot_hessian(self.cache, *outputs)


class REffHam(EffHam):
    """Effective Hamiltonian for restricted wavefunctions."""

    ndm = 1
    deriv_scale = 2.0

    @doc_inherit(EffHam)
    def reset(self, in_dm_alpha):
        self.cache.clear()
        # Take a copy of the input alpha density matrix in the cache.
        dm_alpha = self.cache.load('dm_alpha', alloc=in_dm_alpha.new)[0]
        dm_alpha.assign(in_dm_alpha)

    @doc_inherit(EffHam)
    def reset_delta(self, in_delta_dm_alpha):
        self.cache.clear(tags='d')
        # Take a copy of the input alpha delta density matrix in the cache.
        delta_dm_alpha = self.cache.load('delta_dm_alpha', alloc=in_delta_dm_alpha.new, tags='d')[0]
        delta_dm_alpha.assign(in_delta_dm_alpha)

    @doc_inherit(EffHam)
    def compute_fock(self, fock_alpha):
        EffHam.compute_fock(self, fock_alpha)

    @doc_inherit(EffHam)
    def compute_dot_hessian(self, output_alpha):
        EffHam.compute_dot_hessian(self, output_alpha)


class UEffHam(EffHam):
    """Effective Hamiltonian for unrestricted wavefunctions."""

    ndm = 2

    @doc_inherit(EffHam)
    def reset(self, in_dm_alpha, in_dm_beta):
        self.cache.clear()
        # Take copies of the input alpha and beta density matrices in the cache.
        dm_alpha = self.cache.load('dm_alpha', alloc=in_dm_alpha.new)[0]
        dm_alpha.assign(in_dm_alpha)
        dm_beta = self.cache.load('dm_beta', alloc=in_dm_beta.new)[0]
        dm_beta.assign(in_dm_beta)

    @doc_inherit(EffHam)
    def reset_delta(self, in_delta_dm_alpha, in_delta_dm_beta):
        self.cache.clear(tags='d')
        # Take a copy of the input alpha and beta delta density matrix in the cache.
        delta_dm_alpha = self.cache.load('delta_dm_alpha', alloc=in_delta_dm_alpha.new, tags='d')[0]
        delta_dm_alpha.assign(in_delta_dm_alpha)
        delta_dm_beta = self.cache.load('delta_dm_beta', alloc=in_delta_dm_beta.new, tags='d')[0]
        delta_dm_beta.assign(in_delta_dm_beta)

    @doc_inherit(EffHam)
    def compute_fock(self, fock_alpha, fock_beta):
        EffHam.compute_fock(self, fock_alpha, fock_beta)

    @doc_inherit(EffHam)
    def compute_dot_hessian(self, output_alpha, output_beta):
        EffHam.compute_dot_hessian(self, output_alpha, output_beta)
