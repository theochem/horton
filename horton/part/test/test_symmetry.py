# -*- coding: utf-8 -*-
# Horton is a Density Functional Theory program.
# Copyright (C) 2011-2013 Toon Verstraelen <Toon.Verstraelen@UGent.be>
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


from horton import *
from horton.test.common import get_pentagon_moments


def get_fake_example():
    generators = [
        np.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0]]),
        np.array([[0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 1, 0]]),
    ]
    fracs = np.array([
        [1, 0, 0],
        [0.5, 0.5, 0.5],
    ])
    numbers = np.array([1, 8])
    cell = Cell(None)
    symmetry = Symmetry('fake', generators, fracs, numbers, cell)
    coordinates, numbers, links = symmetry.generate()
    system = System(coordinates, numbers, props={'links': links}, cell=cell)
    assert system.natom == 3
    return system, symmetry


def test_symmetry_scalar():
    system, symmetry = get_fake_example()
    aim_results = {
        'charges': np.array([0.29, 0.31, -0.6]),
        'volumes': np.array([1.2, 1.4, 3.4]),
    }
    sym_results = symmetry_analysis(system, symmetry, aim_results)
    assert len(sym_results) == 2
    stats = sym_results['charges']
    assert abs(stats[:,0] - [0.3, -0.6]).max() < 1e-10
    assert abs(stats[:,1] - [np.std([0.29, 0.31]), 0.0]).max() < 1e-10
    stats = sym_results['volumes']
    assert abs(stats[:,0] - [1.3, 3.4]).max() < 1e-10
    assert abs(stats[:,1] - [np.std([1.2, 1.4]), 0.0]).max() < 1e-10


def test_symmetry_moments():
    system, symmetry = get_fake_example()

    # setup rotated multipole moments
    m0 = get_pentagon_moments()
    m00 = m0.copy()
    m01 = rotate_cartesian_moments(m0, symmetry.generators[1][:,:3])
    m1 = get_pentagon_moments(get_random_rotation())

    # perturb them in a controlled way
    m00[0] += 0.1
    m01[0] -= 0.1
    m00[1] += 0.1
    m01[2] -= 0.1

    # run analysis
    aim_results = {
        'cartesian_moments': np.array([m00, m01, m1]),
    }
    sym_results = symmetry_analysis(system, symmetry, aim_results)

    # check results
    assert len(sym_results) == 1
    stats = sym_results['cartesian_moments']
    assert abs(stats[:,0] - [m0, m1]).max() < 1e-10
    assert abs(stats[1,1]).max() < 1e-10
    assert abs(stats[0,1,:2] - np.std([-0.1, 0.1])).max() < 1e-10
    assert abs(stats[0,1,2:]).max() < 1e-10
