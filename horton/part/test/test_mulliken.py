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

from numpy.testing import assert_almost_equal, assert_equal

from .. mulliken import get_mulliken_operators


def test_mulliken_operators_water_sto3g():
    # data from water_sto3g_hf_g03.fchk
    numbers = np.array([8, 1, 1])
    ncenter = 3
    shell_types = [0, 0, 1, 0, 0]
    shell_maps = [0, 0, 0, 1, 2]
    # density matrix
    dm_full = np.array([[ 2.10503807e+00, -4.39115917e-01,  8.87107806e-02,  7.16046434e-16,
                          6.27282714e-02, -3.27810857e-02, -3.27810800e-02],
                        [-4.39115917e-01,  1.93312061e+00, -5.06512403e-01,  1.03115828e-16,
                         -3.58159960e-01, -1.54725977e-02, -1.54725729e-02],
                        [ 8.87107806e-02, -5.06512403e-01,  1.10208185e+00, -9.19125113e-16,
                          2.69718983e-01,  7.59035717e-02,  6.87997641e-01],
                        [ 7.16046434e-16,  1.03115828e-16, -9.19125113e-16, 2.00000000e+00,
                         -1.38129990e-16, -2.78992150e-16,  1.49460969e-16],
                        [ 6.27282714e-02, -3.58159960e-01,  2.69718983e-01, -1.38129990e-16,
                          9.11364268e-01,  7.02894999e-01, -1.62732921e-01],
                        [-3.27810857e-02, -1.54725977e-02,  7.59035717e-02, -2.78992150e-16,
                          7.02894999e-01,  5.86151074e-01, -1.93687449e-01],
                        [-3.27810800e-02, -1.54725729e-02,  6.87997641e-01, 1.49460969e-16,
                         -1.62732921e-01, -1.93687449e-01,  5.86151122e-01]])
    # overlap matrix
    olp = np.array([[1., 0.23670392, 0., 0., 0., 0.05490733, 0.05490732],
                    [0.23670392, 1., 0., 0., 0., 0.47954331, 0.47954323],
                    [0., 0., 1., 0., 0., 0., 0.37329955],
                    [0., 0., 0., 1., 0., 0., 0.],
                    [0., 0., 0., 0., 1., 0.39594342, -0.13197959],
                    [0.05490733, 0.47954331, 0., 0., 0.39594342, 1., 0.23846113],
                    [0.05490732, 0.47954323, 0.37329955, 0., -0.13197959, 0.23846113, 1.]])
    # compute operators
    operators = get_mulliken_operators(olp, ncenter, shell_types, shell_maps)
    # check operators are symmetric
    for operator in operators:
        assert_equal(operator, operator.T)
    # compute Mulliken charges
    populations = np.array([np.einsum('ab,ba', operator, dm_full) for operator in operators])
    charges = numbers - populations
    # check Mulliken charges
    assert_almost_equal(charges.sum(), 0.0, decimal=6)
    expected = np.array([-3.81897766E-01, 1.90948904E-01, 1.90948861E-01])
    assert_almost_equal(charges, expected, decimal=6)
