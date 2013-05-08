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
'''Symmetry analysis of AIM results'''


import numpy as np

from horton.moments import rotate_cartesian_moments



def symmetry_analysis(system, symmetry, aim_results):
    '''Compute averages and standard deviations on AIM results of equivalent atoms

       **Arguments:**

       system
            The system on which the partitioning was carried out.

       symmetry
            The symmetry descriptor that is used to find the equivalent atoms.

       aim_results
            A dictionary with AIM results. The following fields are supported
            (and the rest gets ignored): 'charges', 'populations',
            'pseudo_populations', 'cartesian_moments', 'radial_moments',
            'volumes', 'volume_ratios', 'c6s'. (See partitioning schemes for
            more details about these fields.)

       **Returns:** a dictionary with the above keys but where each item is an
       array with the corresponding statistical analysis. The dimension of this
       array is one higher:

       * first dimension: number if unique atoms in the primitive unit, instead
         of the number of atoms in the full system
       * second dimension: index can be 0 or 1. 0 refers to the mean over all
         equivalent atoms. 1 refers to the standard deviation
       * remaining dimensions are borrowed from the origin data in aim_results.
    '''
    # define masks to select atoms that match the one atom in the primitive cell.
    links = symmetry.identify(system)
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
            rotated.append(rotate_cartesian_moments(equiv_values[i], generator[:,:3].T))
        return np.mean(rotated, axis=0), np.std(rotated, axis=0)

    # associate analysis functions with certain fields
    cases = [
        (scalar_analysis,
         ['charges', 'populations', 'pseudo_populations',
          'radial_moments', 'volumes', 'volume_ratios', 'c6s']),
        (cartesian_multipole_analysis,
         ['cartesian_moments']),
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
