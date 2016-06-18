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


import numpy as np
from nose.tools import assert_raises
from nose.plugins.attrib import attr

from horton import *  # pylint: disable=wildcard-import,unused-wildcard-import

from horton.test.common import get_random_cell


def check_frac_cart(cell):
    cart1 = np.random.uniform(-20, 20, 3)
    frac1 = cell.to_frac(cart1)
    cart2 = cell.to_cart(frac1)
    frac2 = cell.to_frac(cart2)
    assert abs(cart1 - cart2).max() < 1e-10
    assert abs(frac1 - frac2).max() < 1e-10


def check_g_lincomb_dot_rvecs(cell):
    frac1 = np.random.uniform(-20, 20, 3)
    cart1 = cell.g_lincomb(frac1)
    frac2 = cell.dot_rvecs(cart1)
    cart2 = cell.g_lincomb(frac2)
    assert abs(cart1 - cart2).max() < 1e-10
    assert abs(frac1 - frac2).max() < 1e-10


def test_cell_cubic():
    rvecs = np.array([[9.865, 0.0, 0.0], [0.0, 9.865, 0.0], [0.0, 0.0, 9.865]])*angstrom
    cell = Cell(rvecs)

    # Test attributes
    assert cell.nvec == 3
    assert (cell.rvecs == rvecs).all()
    assert (cell.rspacings == 9.865*angstrom).all()
    assert (cell.gspacings == 1/(9.865*angstrom)).all()
    assert abs(cell.volume - (9.865*angstrom)**3) < 1e-10
    assert abs(np.dot(cell.gvecs, cell.rvecs.transpose()) - np.identity(3)).max() < 1e-5
    assert abs(np.dot(cell.gvecs.transpose(), cell.rvecs) - np.identity(3)).max() < 1e-5
    cell2 = Cell(-cell.rvecs)
    assert abs(cell2.volume - (9.865*angstrom)**3) < 1e-10
    for i in xrange(3):
        assert cell.get_rlength(i) == cell.rlengths[i]
        assert cell.get_glength(i) == cell.glengths[i]
        assert cell.get_rspacing(i) == cell.rspacings[i]
        assert cell.get_gspacing(i) == cell.gspacings[i]
        assert abs(cell.get_rlength(i) - 1.0/cell.get_gspacing(i)) < 1e-10
        assert abs(cell.get_glength(i) - 1.0/cell.get_rspacing(i)) < 1e-10

    # Test methods (1)
    vec1 = np.array([10.0, 0.0, 5.0])*angstrom
    cell.mic(vec1)
    assert abs(vec1 - np.array([0.135, 0.0, -4.865])*angstrom).max() < 1e-10
    cell.add_rvec(vec1, np.array([1,2,3]))
    assert abs(vec1 - np.array([10.0, 19.73, 24.73])*angstrom).max() < 1e-10

    # Test methods (2)
    check_frac_cart(cell)
    check_g_lincomb_dot_rvecs(cell)


def test_cell_triclinic():
    while True:
        rvecs = np.random.uniform(-1, 1, (3, 3))
        if abs(np.linalg.det(rvecs)) > 0.1:
            break
    cell = Cell(rvecs)

    # Test attributes
    assert cell.nvec == 3
    assert (cell.rvecs == rvecs).all()
    assert abs(cell.volume - abs(np.linalg.det(rvecs))) < 1e-10
    assert abs(np.dot(cell.gvecs, cell.rvecs.transpose()) - np.identity(3)).max() < 1e-5
    assert abs(np.dot(cell.gvecs.transpose(), cell.rvecs) - np.identity(3)).max() < 1e-5
    cell2 = Cell(-cell.rvecs)
    assert abs(cell2.volume - abs(np.linalg.det(rvecs))) < 1e-10
    for i in xrange(3):
        assert cell.get_rlength(i) == cell.rlengths[i]
        assert cell.get_glength(i) == cell.glengths[i]
        assert cell.get_rspacing(i) == cell.rspacings[i]
        assert cell.get_gspacing(i) == cell.gspacings[i]
        assert abs(cell.get_rlength(i) - 1.0/cell.get_gspacing(i)) < 1e-10
        assert abs(cell.get_glength(i) - 1.0/cell.get_rspacing(i)) < 1e-10

    # Test methods (1)
    vec1 = np.array([10.0, 0.0, 5.0])*angstrom
    cell.mic(vec1)
    cell.add_rvec(vec1, np.array([1,2,3]))

    # Test methods (2)
    check_frac_cart(cell)
    check_g_lincomb_dot_rvecs(cell)


def test_cell_parallellogram2d():
    rvecs = np.array([[4.922, 0.0, 0.0], [2.462, 4.262, 0.0]])*angstrom
    cell = Cell(rvecs)

    # Test attributes
    assert cell.nvec == 2
    assert (cell.rvecs == rvecs).all()
    assert abs(cell.volume - np.linalg.norm(np.cross(cell.rvecs[0], cell.rvecs[1]))) < 1e-10
    assert abs(np.dot(cell.gvecs, cell.rvecs.transpose()) - np.identity(2)).max() < 1e-5
    for i in xrange(2):
        assert cell.get_rlength(i) == cell.rlengths[i]
        assert cell.get_glength(i) == cell.glengths[i]
        assert cell.get_rspacing(i) == cell.rspacings[i]
        assert cell.get_gspacing(i) == cell.gspacings[i]
        assert abs(cell.get_rlength(i) - 1.0/cell.get_gspacing(i)) < 1e-10
        assert abs(cell.get_glength(i) - 1.0/cell.get_rspacing(i)) < 1e-10
    assert abs(cell.get_rlength(2) - 1.0) < 1e-10
    assert abs(cell.get_glength(2) - 1.0) < 1e-10
    assert abs(cell.get_rspacing(2) - 1.0) < 1e-10
    assert abs(cell.get_gspacing(2) - 1.0) < 1e-10

    # Test methods (1)
    vec1 = np.array([10.0, 0.0, 105.0])*angstrom
    cell.mic(vec1)
    assert abs(vec1 - np.array([0.156, 0.0, 105])*angstrom).max() < 1e-3
    cell.add_rvec(vec1, np.array([1,2]))
    assert abs(vec1 - np.array([10.002, 8.524, 105])*angstrom).max() < 1e-3

    # Test methods (2)
    check_frac_cart(cell)
    check_g_lincomb_dot_rvecs(cell)


def test_cell_1d():
    rvecs = np.array([[5.075, 0.187, 0.055]])*angstrom
    cell = Cell(rvecs)

    # Test attributes
    assert cell.nvec == 1
    assert (cell.rvecs == rvecs).all()
    assert cell.rvecs.shape == (1, 3)
    assert cell.gvecs.shape == (1, 3)
    assert abs(cell.volume - np.linalg.norm(cell.rvecs[0])) < 1e-10
    assert abs(np.dot(cell.gvecs, cell.rvecs.transpose()) - 1) < 1e-5
    assert cell.get_rlength(0) == cell.rlengths[0]
    assert cell.get_glength(0) == cell.glengths[0]
    assert cell.get_rspacing(0) == cell.rspacings[0]
    assert cell.get_gspacing(0) == cell.gspacings[0]
    assert abs(cell.get_rlength(0) - 1.0/cell.get_gspacing(0)) < 1e-10
    assert abs(cell.get_glength(0) - 1.0/cell.get_rspacing(0)) < 1e-10
    for i in xrange(1,3):
        assert abs(cell.get_rlength(i) - 1.0) < 1e-10
        assert abs(cell.get_glength(i) - 1.0) < 1e-10
        assert abs(cell.get_rspacing(i) - 1.0) < 1e-10
        assert abs(cell.get_gspacing(i) - 1.0) < 1e-10
    assert abs(cell.get_rspacing(0)*cell.get_gspacing(0) - 1.0) < 1e-10

    # Test methods (1)
    vec1 = np.array([10.0, 0.0, 105.0])*angstrom
    cell.mic(vec1)
    assert abs(vec1 - np.array([-0.15, -0.374, 104.89])*angstrom).max() < 1e-3
    cell.add_rvec(vec1, np.array([1]))
    assert abs(vec1 - np.array([4.925, -0.187, 104.945])*angstrom).max() < 1e-3

    # Test methods (2)
    check_frac_cart(cell)
    check_g_lincomb_dot_rvecs(cell)


def test_cell_quartz():
    cell = Cell(np.array([[0.0, 0.0, 5.405222], [0.0, 4.913416, 0.0], [-4.255154, 2.456708, 0.0]])*angstrom)

    # Test attributes
    assert cell.rvecs.shape == (3, 3)
    assert cell.gvecs.shape == (3, 3)
    assert abs(cell.volume - abs(np.linalg.det(cell.rvecs))) < 1e-10
    assert abs(np.dot(cell.gvecs, cell.rvecs.transpose()) - np.identity(3)).max() < 1e-5
    for i in xrange(3):
        assert cell.get_rlength(i) == cell.rlengths[i]
        assert cell.get_glength(i) == cell.glengths[i]
        assert cell.get_rspacing(i) == cell.rspacings[i]
        assert cell.get_gspacing(i) == cell.gspacings[i]
        assert abs(cell.get_rlength(i) - 1.0/cell.get_gspacing(i)) < 1e-10
        assert abs(cell.get_glength(i) - 1.0/cell.get_rspacing(i)) < 1e-10

    # Test methods (2)
    check_frac_cart(cell)
    check_g_lincomb_dot_rvecs(cell)

    # Test domain errors
    for i in -1, 4, 245:
        with assert_raises(ValueError):
            cell.get_rspacing(i)
        with assert_raises(ValueError):
            cell.get_gspacing(i)


def test_cell_0d():
    cell = Cell()

    # Test attributes
    assert cell.nvec == 0
    assert cell.rvecs.shape == (0, 3)
    assert cell.gvecs.shape == (0, 3)
    assert cell.rspacings.shape == (0,)
    assert cell.gspacings.shape == (0,)
    for i in xrange(3):
        assert abs(cell.get_rlength(i) - 1.0) < 1e-10
        assert abs(cell.get_glength(i) - 1.0) < 1e-10
        assert abs(cell.get_rspacing(i) - 1.0) < 1e-10
        assert abs(cell.get_gspacing(i) - 1.0) < 1e-10

    # Test methods (1)
    vec1 = np.array([10.0, 0.0, 105.0])*angstrom
    cell.mic(vec1)
    assert abs(vec1 - np.array([10.0, 0.0, 105.0])*angstrom).max() < 1e-3
    cell.add_rvec(vec1, np.array([], dtype=int))
    assert abs(vec1 - np.array([10.0, 0.0, 105.0])*angstrom).max() < 1e-3

    # Test methods (2)
    check_frac_cart(cell)
    check_g_lincomb_dot_rvecs(cell)


def setup_ranges_rcut(nvec):
    a = 10**np.random.uniform(-1, 1)
    cell = get_random_cell(a, nvec)

    origin = np.random.uniform(-3*a, 3*a, 3)
    center = np.random.uniform(-3*a, 3*a, 3)
    rcut = np.random.uniform(0.2*a, 5*a)

    ranges_begin, ranges_end = cell.get_ranges_rcut(origin-center, rcut)
    ranges_low = ranges_begin-2
    ranges_high = ranges_end+2

    return cell, origin, center, rcut, ranges_begin, ranges_end, ranges_low, ranges_high


@attr('slow')
def test_ranges_rcut_3d():
    counter = 0
    while True:
        cell, origin, center, rcut, ranges_begin, ranges_end, ranges_low, ranges_high = setup_ranges_rcut(3)
        npoint = np.product(ranges_high-ranges_low)
        if npoint > 10000:
            continue

        # test if distances to points outside the ranges are always outside rcut
        for i0 in xrange(ranges_low[0], ranges_high[0]):
            for i1 in xrange(ranges_low[1], ranges_high[1]):
                for i2 in xrange(ranges_low[2], ranges_high[2]):
                    if (i0 >= ranges_begin[0] and i0 < ranges_end[0] and
                        i1 >= ranges_begin[1] and i1 < ranges_end[1] and
                        i2 >= ranges_begin[2] and i2 < ranges_end[2]):
                        continue
                    tmp = origin.copy()
                    cell.add_rvec(tmp, np.array([i0, i1, i2]))
                    assert np.linalg.norm(tmp - center) > rcut

        counter += 1
        if counter >= 20:
            break


def test_ranges_rcut_2d():
    counter = 0
    while True:
        cell, origin, center, rcut, ranges_begin, ranges_end, ranges_low, ranges_high = setup_ranges_rcut(2)

        # test if distances to points outside the ranges are always outside rcut
        for i0 in xrange(ranges_low[0], ranges_high[0]):
            for i1 in xrange(ranges_low[1], ranges_high[1]):
                if (i0 >= ranges_begin[0] and i0 < ranges_end[0] and
                    i1 >= ranges_begin[1] and i1 < ranges_end[1]):
                    continue
                tmp = origin.copy()
                cell.add_rvec(tmp, np.array([i0, i1]))
                assert np.linalg.norm(tmp - center) > rcut

        counter += 1
        if counter >= 50:
            break


def test_ranges_rcut_1d():
    counter = 0
    while True:
        cell, origin, center, rcut, ranges_begin, ranges_end, ranges_low, ranges_high = setup_ranges_rcut(1)

        # test if distances to points outside the ranges are always outside rcut
        for i0 in xrange(ranges_low[0], ranges_high[0]):
            if i0 >= ranges_begin[0] and i0 < ranges_end[0]:
                continue
            tmp = origin.copy()
            cell.add_rvec(tmp, np.array([i0]))
            assert np.linalg.norm(tmp - center) > rcut

        counter += 1
        if counter >= 200:
            break


def test_ranges_rcut_0d():
    cell = Cell(None)
    a = 0.25
    origin = np.random.uniform(-3*a, 3*a, 3)
    center = np.random.uniform(-3*a, 3*a, 3)
    rcut = np.random.uniform(0.2*a, 5*a)

    ranges_begin, ranges_end = cell.get_ranges_rcut(origin-center, rcut)
    assert ranges_begin.size == 0
    assert ranges_begin.shape == (0,)
    assert ranges_end.size == 0
    assert ranges_end.shape == (0,)


def test_smart_wrap():
    assert smart_wrap(3, 5, 0) == 3
    assert smart_wrap(7, 5, 0) == -1
    assert smart_wrap(-2, 5, 0) == -1
    assert smart_wrap(3, 5, 1) == 3
    assert smart_wrap(7, 5, 1) == 2
    assert smart_wrap(-2, 5, 1) == 3


def test_from_parameters0():
    assert Cell.from_parameters([], []).nvec == 0


def check_from_parameters(cell0):
    lengths0, angles0 = cell0.parameters
    cell1 = Cell.from_parameters(lengths0, angles0)
    lengths1, angles1 = cell1.parameters
    assert lengths0.shape == lengths1.shape
    assert abs(lengths0 - lengths1).max() < 1e-10
    assert angles0.shape == angles1.shape
    if angles0.size > 0:
        assert abs(angles0 - angles1).max() < 1e-10

def test_from_parameters1():
    for i in xrange(10):
        cell0 = get_random_cell(1.0, 1)
        check_from_parameters(cell0)


def test_from_parameters2():
    for i in xrange(10):
        cell0 = get_random_cell(1.0, 2)
        check_from_parameters(cell0)


def test_from_parameters3():
    for i in xrange(10):
        cell0 = get_random_cell(1.0, 3)
        check_from_parameters(cell0)


def test_no_initvoid():
    with assert_raises(TypeError):
        cell = Cell(initvoid=True)
