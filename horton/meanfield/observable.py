# -*- coding: utf-8 -*-
# HORTON: Helpful Open-source Research TOol for N-fermion systems.
# Copyright (C) 2011-2022 The HORTON Development Team
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
"""Base classes for energy terms and other observables of the wavefunction"""


import numpy as np

from horton.utils import doc_inherit


__all__ = [
    'compute_dm_full',
    'Observable',
    'RTwoIndexTerm', 'UTwoIndexTerm',
    'RDirectTerm', 'UDirectTerm',
    'RExchangeTerm', 'UExchangeTerm',
]


def compute_dm_full(cache, prefix='', tags=None):
    """Add the spin-summed density matrix to the cache unless it is already present.

    Parameters
    ----------
    cache : Cache
        Used to store intermediate results that can be reused or inspected later.
    prefix : str
        A prefix, in case one wants to load a special density matrix and store the results
        in a special place. The typical use case is the `delta_` prefix which is used to
        process a change in density.
    tags : str
        The tags to use for the cache when allocating the arrays. In case of changes in
        density matrices, this argument is typically equal to 'd'.
    """
    dm_alpha = cache['%sdm_alpha' % prefix]
    dm_beta = cache['%sdm_beta' % prefix]
    dm_full, new = cache.load('%sdm_full' % prefix, alloc=dm_alpha.shape, tags=tags)
    if new:
        dm_full[:] = dm_alpha
        dm_full += dm_beta
    return dm_full


class Observable(object):
    """Base class for contribution to EffHam classes.

    These are usually energy expressions (as function of one or more density matrices).
    One may also use this for other observables, e.g. to construct a Lagrangian instead of
    a regular effective Hamiltonian.
    """

    def __init__(self, label):
        """Initialize an Observable instance.

        Parameters
        ----------
        label : str
            A short string to identify the observable.
        """
        self.label = label

    def compute_energy(self, cache):
        """Compute the expectation value of the observable.

        Parameters
        ----------
        cache : Cache
            Used to store intermediate results that can be reused or inspected later.
        """
        raise NotImplementedError

    def add_fock(self, cache, *focks):
        """Add contributions to the Fock matrices.

        Parameters
        ----------
        cache : Cache
            Used to store intermediate results that can be reused or inspected later.
        fock1, fock2, ... : np.ndarray, shape=(nbasis, nbasis)
            A list of output Fock operators. The caller is responsible for setting these
            operators initially to zero (if desired).
        """
        raise NotImplementedError

    def add_dot_hessian(self, cache, *outputs):
        """Add contribution to the dot product of the Hessian with the trial delta DM.

        The Hessian in this method is the second derivative of the observable towards the
        matrix elements of the density matrix or matrices. The ``dms`` and ``delta dms``
        are stored in the cache.

        Parameters
        ----------
        cache
            A cache object used to store intermediate results that can be reused or
            inspected later. All results specific for the detla_dm are stored with the tag
            'd' in the cache object.
        output1, output2, ... : np.ndarray, shape=(nbasis, nbasis)
            Output objects, to which contributions from the dot product of the Hessian
            with the delta density matrices is added.
        """
        raise NotImplementedError


class RTwoIndexTerm(Observable):
    """Observable linear in the density matrix (restricted)."""

    def __init__(self, op_alpha, label):
        """Initialize a RTwoIndexTerm instance.

        Parameters
        ----------
        op_alpha : np.ndarray, shape=(nbasis, nbasis)
            Expansion of one-body operator in basis of alpha orbitals. Same is used for
            beta.
        label : str
            A short string to identify the observable.
        """
        self.op_alpha = op_alpha
        Observable.__init__(self, label)

    @doc_inherit(Observable)
    def compute_energy(self, cache):
        return 2 * np.einsum('ab,ba', self.op_alpha, cache['dm_alpha'])

    @doc_inherit(Observable)
    def add_fock(self, cache, fock_alpha):
        fock_alpha += self.op_alpha

    @doc_inherit(Observable)
    def add_dot_hessian(self, cache, output_alpha):
        pass


class UTwoIndexTerm(Observable):
    """Observable linear in the density matrix (unrestricted)."""

    def __init__(self, op_alpha, label, op_beta=None):
        """Initialize a RTwoIndexTerm instance.

        Parameters
        ----------
        op_alpha : np.ndarray, shape=(nbasis, nbasis)
            Expansion of one-body operator in basis of alpha orbitals.
        label : str
            A short string to identify the observable.
        op_beta : np.ndarray, shape=(nbasis, nbasis)
            Expansion of one-body operator in basis of beta orbitals. When not given,
            op_alpha is used.
        """
        self.op_alpha = op_alpha
        self.op_beta = op_alpha if op_beta is None else op_beta
        Observable.__init__(self, label)

    @doc_inherit(Observable)
    def compute_energy(self, cache):
        if self.op_alpha is self.op_beta:
            # when both operators are references to the same object, take a
            # shortcut
            compute_dm_full(cache)
            return np.einsum('ab,ba', self.op_alpha, cache['dm_full'])
        else:
            # If the operator is different for different spins, do the normal
            # thing.
            return np.einsum('ab,ba', self.op_alpha, cache['dm_alpha']) + \
                   np.einsum('ab,ba', self.op_beta, cache['dm_beta'])

    @doc_inherit(Observable)
    def add_fock(self, cache, fock_alpha, fock_beta):
        fock_alpha += self.op_alpha
        fock_beta += self.op_beta

    @doc_inherit(Observable)
    def add_dot_hessian(self, cache, output_alpha, output_beta):
        pass


def contract_direct(op, dm):
    """Perform an direct-type contraction with a four-index operator.

    Parameters
    ----------
    op : np.ndarray, shape=(nvec, nbasis, nbasis) or
        The four-index operator or its Cholesky decomposition
    dm : np.ndarray, shape=(nbasis, nbasis)
        The density matrix
    """
    if op.ndim == 3:
        # Cholesky decomposition
        tmp = np.tensordot(op, dm, axes=([(1, 2), (1, 0)]))
        return np.tensordot(op, tmp, [0, 0])
    elif op.ndim == 4:
        # Normal case
        return np.einsum('abcd,bd->ac', op, dm)
    else:
        raise NotImplementedError


class RDirectTerm(Observable):
    """Direct term of the expectation value of a two-body operator (restricted)."""

    def __init__(self, op_alpha, label):
        """Initialize a RDirectTerm instance.

        Parameters
        ----------
        op_alpha
            Expansion of two-body operator in basis of alpha orbitals. Same is used for
            beta orbitals. Also a Cholesky decomposition of the operator is supported.
        label : str
            A short string to identify the observable.
        """
        self.op_alpha = op_alpha
        Observable.__init__(self, label)

    def _update_direct(self, cache):
        """Recompute the direct operator if it has become invalid.

        Parameters
        ----------
        cache : Cache
            Used to store intermediate results that can be reused or inspected later.
        """
        dm_alpha = cache['dm_alpha']
        direct, new = cache.load('op_%s_alpha' % self.label, alloc=dm_alpha.shape)
        if new:
            direct[:] = contract_direct(self.op_alpha, dm_alpha)
            direct *= 2  # contribution from beta electrons is identical

    @doc_inherit(Observable)
    def compute_energy(self, cache):
        self._update_direct(cache)
        direct = cache.load('op_%s_alpha' % self.label)
        return np.einsum('ab,ba', direct, cache['dm_alpha'])

    @doc_inherit(Observable)
    def add_fock(self, cache, fock_alpha):
        self._update_direct(cache)
        direct = cache.load('op_%s_alpha' % self.label)
        fock_alpha += direct

    @doc_inherit(Observable)
    def add_dot_hessian(self, cache, output_alpha):
        delta_dm_alpha = cache.load('delta_dm_alpha')
        output_alpha += contract_direct(self.op_alpha, delta_dm_alpha)


class UDirectTerm(Observable):
    """Direct term of the expectation value of a two-body operator (unrestricted)."""

    def __init__(self, op_alpha, label, op_beta=None):
        """Initialize a UDirectTerm instance.

        Parameters
        ----------
        op_alpha
            Expansion of two-body operator in basis of alpha orbitals. Also a Cholesky
            decomposition of the operator is supported.
        label : str
            A short string to identify the observable.
        op_beta
            Expansion of two-body operator in basis of beta orbitals. When not given,
            op_alpha is used. Also a Cholesky decomposition of the operator is supported.
        """
        self.op_alpha = op_alpha
        self.op_beta = op_alpha if op_beta is None else op_beta
        Observable.__init__(self, label)

    def _update_direct(self, cache):
        """Recompute the direct operator(s) if it/they has/have become invalid.

        Parameters
        ----------
        cache : Cache
            Used to store intermediate results that can be reused or inspected later.
        """
        if self.op_alpha is self.op_beta:
            # This branch is nearly always going to be followed in practice.
            dm_full = compute_dm_full(cache)
            direct, new = cache.load('op_%s' % self.label, alloc=dm_full.shape)
            if new:
                direct[:] = contract_direct(self.op_alpha, dm_full)
        else:
            # This is probably never going to happen. In case it does, please
            # add the proper code here.
            raise NotImplementedError

    @doc_inherit(Observable)
    def compute_energy(self, cache):
        self._update_direct(cache)
        if self.op_alpha is self.op_beta:
            # This branch is nearly always going to be followed in practice.
            direct = cache['op_%s' % self.label]
            dm_full = cache['dm_full']
            return 0.5 * np.einsum('ab,ab', direct, dm_full)
        else:
            # This is probably never going to happen. In case it does, please
            # add the proper code here.
            raise NotImplementedError

    @doc_inherit(Observable)
    def add_fock(self, cache, fock_alpha, fock_beta):
        self._update_direct(cache)
        if self.op_alpha is self.op_beta:
            # This branch is nearly always going to be followed in practice.
            direct = cache['op_%s' % self.label]
            fock_alpha += direct
            fock_beta += direct
        else:
            # This is probably never going to happen. In case it does, please
            # add the proper code here.
            raise NotImplementedError

    @doc_inherit(Observable)
    def add_dot_hessian(self, cache, output_alpha, output_beta):
        if self.op_alpha is self.op_beta:
            delta_dm_full = compute_dm_full(cache, prefix='delta_', tags='d')
            output_alpha += contract_direct(self.op_alpha, delta_dm_full)
            output_beta += contract_direct(self.op_alpha, delta_dm_full)
        else:
            # This is probably never going to happen. In case it does, please
            # add the proper code here.
            raise NotImplementedError


def contract_exchange(op, dm):
    """Perform an exchange-type contraction with a four-index operator.

    Parameters
    ----------
    op : np.ndarray, shape=(nvec, nbasis, nbasis) or
        The four-index operator or its Cholesky decomposition
    dm : np.ndarray, shape=(nbasis, nbasis)
        The density matrix
    """
    if op.ndim == 3:
        # Cholesky decomposition
        tmp = np.tensordot(op, dm, axes=([1, 1]))
        return np.tensordot(op, tmp, ([0, 2], [0, 2]))
    elif op.ndim == 4:
        return np.einsum('abcd,cb->ad', op, dm)
    else:
        raise NotImplementedError


class RExchangeTerm(Observable):
    """Exchange term of the expectation value of a two-body operator (restricted)."""

    def __init__(self, op_alpha, label, fraction=1.0):
        """Initialize a RExchangeTerm instance.

        Parameters
        ----------
        op_alpha
            Expansion of two-body operator in basis of alpha orbitals. Same is used for
            beta orbitals. Also a Cholesky decomposition of the operator is supported.
        fraction : float
            Amount of exchange to be included (1.0 corresponds to 100%).
        label : str
            A short string to identify the observable.
        """
        self.op_alpha = op_alpha
        self.fraction = fraction
        Observable.__init__(self, label)

    def _update_exchange(self, cache):
        """Recompute the Exchange operator if invalid.

        Parameters
        ----------
        cache : Cache
            Used to store intermediate results that can be reused or inspected later.
        """
        dm_alpha = cache['dm_alpha']
        exchange_alpha, new = cache.load('op_%s_alpha' % self.label,
                                         alloc=dm_alpha.shape)
        if new:
            exchange_alpha[:] = contract_exchange(self.op_alpha, dm_alpha)

    @doc_inherit(Observable)
    def compute_energy(self, cache):
        self._update_exchange(cache)
        exchange_alpha = cache['op_%s_alpha' % self.label]
        dm_alpha = cache['dm_alpha']
        return -self.fraction * np.einsum('ab,ba', exchange_alpha, dm_alpha)

    @doc_inherit(Observable)
    def add_fock(self, cache, fock_alpha):
        self._update_exchange(cache)
        exchange_alpha = cache['op_%s_alpha' % self.label]
        fock_alpha -= self.fraction*exchange_alpha

    @doc_inherit(Observable)
    def add_dot_hessian(self, cache, output_alpha):
        delta_dm_alpha = cache.load('delta_dm_alpha')
        output_alpha -= (0.5*self.fraction)*contract_exchange(self.op_alpha, delta_dm_alpha)


class UExchangeTerm(Observable):
    """Exchange term of the expectation value of a two-body operator (unrestricted)."""

    def __init__(self, op_alpha, label, fraction=1.0, op_beta=None):
        """Initialize a UExchangeTerm instance.

        Parameters
        ----------
        op_alpha
            Expansion of two-body operator in basis of alpha orbitals. Also a Cholesky
            decomposition of the operator is supported.
        label : str
            A short string to identify the observable.
        fraction : float
            Amount of exchange to be included (1.0 corresponds to 100%).
        op_beta
            Expansion of two-body operator in basis of beta orbitals. When not given,
            op_alpha is used. Also a Cholesky decomposition of the operator is supported.
        """
        self.op_alpha = op_alpha
        self.op_beta = op_alpha if op_beta is None else op_beta
        self.fraction = fraction
        Observable.__init__(self, label)

    def _update_exchange(self, cache):
        """Recompute the Exchange operator(s) if invalid.

        Parameters
        ----------
        cache : Cache
            Used to store intermediate results that can be reused or inspected later.
        """
        # alpha
        dm_alpha = cache['dm_alpha']
        exchange_alpha, new = cache.load('op_%s_alpha' % self.label,
                                         alloc=dm_alpha.shape)
        if new:
            exchange_alpha[:] = contract_exchange(self.op_alpha, dm_alpha)
        # beta
        dm_beta = cache['dm_beta']
        exchange_beta, new = cache.load('op_%s_beta' % self.label,
                                         alloc=dm_beta.shape)
        if new:
            exchange_beta[:] = contract_exchange(self.op_beta, dm_beta)

    @doc_inherit(Observable)
    def compute_energy(self, cache):
        self._update_exchange(cache)
        exchange_alpha = cache['op_%s_alpha' % self.label]
        exchange_beta = cache['op_%s_beta' % self.label]
        dm_alpha = cache['dm_alpha']
        dm_beta = cache['dm_beta']
        return (-0.5*self.fraction)*np.einsum('ab,ba', exchange_alpha, dm_alpha) + \
               (-0.5*self.fraction)*np.einsum('ab,ba', exchange_beta, dm_beta)

    @doc_inherit(Observable)
    def add_fock(self, cache, fock_alpha, fock_beta):
        self._update_exchange(cache)
        exchange_alpha = cache['op_%s_alpha' % self.label]
        fock_alpha -= self.fraction*exchange_alpha
        exchange_beta = cache['op_%s_beta' % self.label]
        fock_beta -= self.fraction*exchange_beta

    @doc_inherit(Observable)
    def add_dot_hessian(self, cache, output_alpha, output_beta):
        delta_dm_alpha = cache.load('delta_dm_alpha')
        output_alpha -= self.fraction*contract_exchange(self.op_alpha, delta_dm_alpha)
        delta_dm_beta = cache.load('delta_dm_beta')
        output_beta -= self.fraction*contract_exchange(self.op_beta, delta_dm_beta)
