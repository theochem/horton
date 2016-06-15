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
'''Utility functions for the ``horton-wpart.py`` script'''

import numpy as np

from horton.log import log
from horton.moments import get_npure_cumul
from horton.cext import fill_pure_polynomials
from horton.meanfield.response import compute_noninteracting_response
from horton.meanfield.bond_order import compute_bond_orders_cs, \
    compute_bond_orders_os

import horton.part
from horton.part.base import WPart

__all__ = ['wpart_schemes', 'wpart_slow_analysis']


def get_wpart_schemes():
    """Return a dictionary with all wpart schemes"""
    wpart_schemes = {}
    for o in vars(horton.part).itervalues():
        if isinstance(o, type) and issubclass(o, WPart) and o.name is not None:
            wpart_schemes[o.name] = o
    return wpart_schemes


wpart_schemes = get_wpart_schemes()


def wpart_slow_analysis(wpart, mol):
    """An additional and optional analysis for horton-wpart.py

       This analysis is currently not included in horton/part because it would
       imply additional changes in the API.

       **Arguments:**

       wpart
            An instance of a WPart class

       mol
            An instance of IOData. This instance must at least contain an
            obasis, exp_alpha and exp_beta (in case of unrestricted spin) object.
    """

    # A) Compute AIM overlap operators
    # --------------------------------
    # These are not included in the output by default because they occupy too
    # much space.
    wpart.do_partitioning()
    npure = get_npure_cumul(wpart.lmax)
    for index in xrange(wpart.natom):
        if log.do_medium:
            log('Computing overlap matrices for atom %i.' % index)

        # Prepare solid harmonics on grids.
        grid = wpart.get_grid(index)
        if wpart.lmax > 0:
            work = np.zeros((grid.size, npure - 1), float)
            work[:, 0] = grid.points[:, 2] - wpart.coordinates[index, 2]
            work[:, 1] = grid.points[:, 0] - wpart.coordinates[index, 0]
            work[:, 2] = grid.points[:, 1] - wpart.coordinates[index, 1]
            if wpart.lmax > 1:
                fill_pure_polynomials(work, wpart.lmax)
        at_weights = wpart.cache.load('at_weights', index)

        # Convert the weight functions to AIM overlap operators.
        counter = 0
        overlap_operators = {}
        for l in xrange(wpart.lmax + 1):
            for m in xrange(-l, l + 1):
                op = mol.lf.create_two_index()
                if counter > 0:
                    tmp = at_weights * work[:, counter - 1]
                else:
                    tmp = at_weights
                mol.obasis.compute_grid_density_fock(grid.points, grid.weights,
                                                     tmp, op)
                overlap_operators['olp_%05i' % counter] = op
                counter += 1

        wpart.cache.dump(('overlap_operators', index), overlap_operators)

    # Correct the s-type overlap operators such that the sum is exactly
    # equal to the total overlap.
    error_overlap = mol.lf.create_two_index()
    for index in xrange(wpart.natom):
        atom_overlap = wpart.cache.load('overlap_operators', index)['olp_00000']
        error_overlap.iadd(atom_overlap)
    error_overlap.iadd(mol.obasis.compute_overlap(mol.lf), -1)
    error_overlap.iscale(1.0 / wpart.natom)
    for index in xrange(wpart.natom):
        atom_overlap = wpart.cache.load('overlap_operators', index)['olp_00000']
        atom_overlap.iadd(error_overlap, -1)

    # A') Construct the list of operators with a logical ordering
    #   * outer loop: s, pz, px, py, ...
    #   * inner loop: atoms
    operators = []
    for ipure in xrange(npure):
        for iatom in xrange(wpart.natom):
            operators.append(wpart.cache.load('overlap_operators', iatom)['olp_%05i' % ipure])

    # B) Compute Wiberg bond orders from the first-order density matrix
    # -----------------------------------------------------------------
    # This is inherently limited because the second-order density matrix is
    # required to compute a proper density matrix. For a pure single-reference
    # method like HF, this is fine. In case of DFT, the interpretation of the
    # bond-orders is dicy. In case of Post-HF, these bond orders are plain
    # wrong.
    if log.do_medium:
        log('Computing bond orders.')
    dm_full = mol.get_dm_full()
    dm_spin = mol.get_dm_spin()
    if dm_spin is None:
        # closed-shell case
        dm_alpha = dm_full.copy()
        dm_alpha.iscale(0.5)
        bond_orders, valences, free_valences = compute_bond_orders_cs(dm_alpha,
                                                                      operators)
    else:
        # open-shell case
        dm_alpha = dm_full.copy()
        dm_alpha.iadd(dm_spin)
        dm_alpha.iscale(0.5)
        dm_beta = dm_full.copy()
        dm_beta.iadd(dm_spin, -1.0)
        dm_beta.iscale(0.5)
        bond_orders, valences, free_valences = compute_bond_orders_os(dm_alpha,
                                                                      dm_beta,
                                                                      operators)
    wpart.cache.dump('bond_orders', bond_orders, tags='o')
    wpart.cache.dump('valences', valences, tags='o')
    wpart.cache.dump('free_valences', free_valences, tags='o')

    # C) Non-interacting response
    # ---------------------------
    # The current implementation is only sensible for a single-reference method
    # as it is based on HF/KS orbitals. The results are meaningful for both
    # HF and KS orbitals.
    if log.do_medium:
        log('Computing Xs response.')
    if dm_spin is None:
        if not hasattr(mol, 'exp_alpha'):
            return
        xs_response = 2 * compute_noninteracting_response(mol.exp_alpha, operators)
    else:
        if not (hasattr(mol, 'exp_alpha') and hasattr(mol, 'exp_beta')):
            return
        xs_response = compute_noninteracting_response(mol.exp_alpha, operators)
        xs_response += compute_noninteracting_response(mol.exp_beta, operators)
    wpart.cache.dump('noninteracting_response', xs_response, tags='o')
