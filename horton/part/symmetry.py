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
'''Symmetry analysis of atoms-in-molecules results'''


import numpy as np

from horton.moments import rotate_cartesian_moments_all


__all__ = ['symmetry_analysis']


def symmetry_analysis(coordinates, cell, symmetry, aim_results):
    '''Compute averages and standard deviations on AIM results of equivalent atoms

       **Arguments:**

       coordinates
            An (N, 3) array of atomic coordinates that adhere (with some minor
            deviation) to this symmetry.

       cell
            A Cell instance describing the periodic boundary conditions

       symmetry
            The symmetry descriptor that is used to find the equivalent atoms.

       aim_results
            A dictionary with AIM results. The following fields are supported
            (and the rest gets ignored): 'charges', 'populations',
            'pseudo_populations', 'cartesian_multipoles', 'radial_moments',
            'volumes', 'volume_ratios', 'c6s'. (See partitioning schemes for
            more details about these fields.)

       **Returns:** a dictionary with the above keys but where each item is an
       array with the corresponding statistical analysis. The dimension of this
       array is one higher:

       * first dimension: number if unique atoms in the primitive unit, instead
         of the number of atoms in the full molecule
       * second dimension: index can be 0 or 1. 0 refers to the mean over all
         equivalent atoms. 1 refers to the standard deviation
       * remaining dimensions are borrowed from the origin data in aim_results.
    '''
    # define masks to select atoms that match the one atom in the primitive cell.
    links = symmetry.identify(coordinates, cell)
    masks = []
    generators = []
    for iprim in xrange(symmetry.natom):
        mask = links[:,0] == iprim
        masks.append(mask)
        generators.append([symmetry.generators[j] for j in links[:,1][mask]])

    # define analysis helpers for scalar and cartesian multipole quantities:
    def scalar_analysis(equiv_values, equiv_generators):
        return equiv_values.mean(axis=0), equiv_values.std(axis=0)

    def cartesian_multipole_analysis(equiv_values, equiv_generators):
        rotated = []
        for i, generator in enumerate(equiv_generators):
            # The multipole moments are rotated back to the primitive unit,
            # so we need to take the transpose.
            rotated.append(rotate_cartesian_moments_all(generator[:,:3].T, equiv_values[i]))
        return np.mean(rotated, axis=0), np.std(rotated, axis=0)

    # associate analysis functions with certain fields
    cases = [
        (scalar_analysis,
         ['charges', 'populations', 'pseudo_populations', 'spin_charges',
          'radial_moments', 'volumes', 'volume_ratios', 'c6s']),
        (cartesian_multipole_analysis,
         ['cartesian_multipoles']),
    ]

    # process everything
    result = {}
    for analysis_fn, fields in cases:
        for key in fields:
            sys_values = aim_results.get(key)
            if sys_values is None:
                continue

            # loop over the atoms in the primitive unit.
            stats = np.zeros((symmetry.natom, 2,) + sys_values.shape[1:], float)
            for iprim in xrange(symmetry.natom):
                equiv_values = sys_values[masks[iprim]]
                stats[iprim] = analysis_fn(equiv_values, generators[iprim])
            result[key] = stats

    return result
