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
r'''Generic implementation of bond orders for mean-field wavefunctions

   In the context bond orders and self-electron delocatization indices (SEDI's)
   are always one an the same thing. For two AIM overlap perators, :math:`S_A`
   and :math:`S_B`, the bond order is defined as:

   .. math::
       \text{BO}_{AB} = 2 \mathrm{Tr}\left[
            (D^{\alpha} S^A)(D^{\alpha} S^B) + (D^{\beta} S^A)(D^{\beta} S^B)
       \right]

   where :math:`D^{\alpha}` and :math:`D^{\beta}` are the density matrices of the
   :math:`\alpha` and :math:`\beta` electrons, respectively. A related quantity
   is the valence of an atom. It is defined as:

   .. math::
       V_A = 2 N_A - \mathrm{Tr}\left[ (D S^A)(D S^A) \right]

   where :math:`D` is the sum of the :math:`\alpha` and :math:`\beta` electron
   density matrices, and :math:`N_A` is the population of atom :math:`A`. The
   free valence is defined as:

   .. math::
       F_A = V_A - \sum_{B \neq A} \text{BO}_{AB}
'''


import numpy as np


__all__ = ['compute_bond_orders_cs', 'compute_bond_orders_os']


def compute_bond_orders_cs(dm_alpha, operators):
    '''Compute bond orders, valences and free valences (closed-shell case)

       **Arguments:**

       dm_alpha
            The density matrix of the alpha electrons

       operators
            A list of one-body operators.

       **Returns:**

       bond_orders
            A symmetric N by N matrix with bond orders.

       valences
            A vector with atomic valences

       free_valences
            A vector with atomic free valences
    '''
    bond_orders, populations = _compute_bond_orders_low(dm_alpha, operators)
    bond_orders *= 2
    populations *= 2
    valences = 4*_compute_valences_low(dm_alpha, populations/4, operators)
    free_valences = valences - (bond_orders.sum(axis=0) - np.diag(bond_orders))
    return bond_orders, valences, free_valences


def compute_bond_orders_os(dm_alpha, dm_beta, operators):
    '''Compute bond orders, valences and free valences (open-shell case)

       **Arguments:**

       dm_alpha
            The density matrix of the alpha electrons

       operators
            A list of one-body operators.

       **Returns:**

       bond_orders
            A symmetric N by N matrix with bond orders.

       valences
            A vector with atomic valences

       free_valences
            A vector with atomic free valences
    '''
    bond_orders_a, populations_a = _compute_bond_orders_low(dm_alpha, operators)
    bond_orders_b, populations_b = _compute_bond_orders_low(dm_beta, operators)
    bond_orders = bond_orders_a + bond_orders_b
    populations = populations_a + populations_b
    dm_full = dm_alpha.copy()
    dm_full.iadd(dm_beta)
    valences = _compute_valences_low(dm_full, populations, operators)
    free_valences = valences - (bond_orders.sum(axis=0) - np.diag(bond_orders))
    return bond_orders, valences, free_valences


def _compute_bond_orders_low(dm, operators):
    '''Compute bond orders and populations

       **Arguments:**

       dm
            An instance of OneBody with a single-electron density matrix

       operators
            A list of one-body operators.

       **Returns:**

       bond_orders
            A symmetric N by N matrix with bond orders.

       populations
            A vector with atomic populations
    '''
    # Run some initial tests
    for op in operators:
        if op.nbasis != dm.nbasis:
            raise TypeError('Mismatch detected between nbasis in op and dm.')

    n = len(operators)
    bond_orders = np.zeros((n, n), float)
    populations = np.zeros(n, float)

    precomputed = []
    for i0 in xrange(n):
        # compute population
        populations[i0] = operators[i0].contract_two('ab,ab', dm)
        # precompute a dot product
        tmp = dm.copy()
        tmp.idot(operators[i0])
        precomputed.append(tmp)
        for i1 in xrange(i0+1):
            # compute bond orders
            bond_orders[i0, i1] = 2*precomputed[i0].contract_two('ab,ba', precomputed[i1])
            bond_orders[i1, i0] = bond_orders[i0, i1]

    return bond_orders, populations


def _compute_valences_low(dm, populations, operators):
    '''Computes the valences

       **Arguments:**

       dm
            An instance of OneBody with the sum of the alpha and beta density
            matrices.

       populations
            The expectation values of the operators for the full density matrix.

       operators
            A list of one-body operators.
    '''
    # Run some initial tests
    for op in operators:
        if op.nbasis != dm.nbasis:
            raise TypeError('Mismatch detected between nbasis in op and dm.')

    n = len(operators)
    valences = np.zeros(n, float)
    for i in xrange(n):
        # valence for atom i
        tmp = dm.copy()
        tmp.idot(operators[i])
        valences[i] = 2*populations[i] - tmp.contract_two('ab,ba', tmp)
    return valences
