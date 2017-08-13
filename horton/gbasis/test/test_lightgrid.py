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


import numpy as np

from .lightgrid import generate_molecular_grid


def test_gaussians_triatomic():
    # An L-shaped CH2 molecule.
    coordinates = np.array([[1.6, 0.0, 0.0], [0.0, 1.6, 0.0], [0.0, 0.0, 0.0]])
    numbers = np.array([1.0, 1.0, 6.0])
    points, weights = generate_molecular_grid(numbers, coordinates, 10000)

    # Put normalized Slaters on the nuclei
    rho = 0.0
    alphas = [2.0, 2.0, 1.8]
    for iatom in range(3):
        deltas = points - coordinates[iatom]
        distances = np.sqrt(np.einsum('ij,ij->i', deltas, deltas))
        rho += (alphas[iatom]**3/(8*np.pi))*np.exp(-distances*alphas[iatom])

    # The test
    np.testing.assert_almost_equal(np.dot(rho, weights), 3.0, decimal=2)
