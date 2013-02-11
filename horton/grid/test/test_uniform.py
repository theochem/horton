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


import numpy as np

from horton import *
from horton.grid.test.common import *


def test_uig_gauss():
    # grid parameters
    spacing = 0.1
    naxis = 81

    # grid setup
    offset = 0.5*spacing*(naxis-1)
    origin = np.zeros(3, float)-offset
    cell = Cell(np.identity(3)*spacing)
    shape = np.array([naxis, naxis, naxis])
    uig = UniformIntGrid(origin, cell, shape)

    # fill a 3D grid with a gaussian function
    rvecs = cell.rvecs
    x, y, z = np.indices(shape)
    dsq = (x*spacing-offset)**2 + (y*spacing-offset)**2 + (z*spacing-offset)**2
    data = np.exp(-0.5*dsq)/(2*np.pi)**(1.5)

    # compute the integral and compare with analytical result
    num_result = uig.integrate(data)
    assert abs(num_result - 1.0) < 1e-3


def test_index_wrap():
    assert index_wrap(2, 3) == 2
    assert index_wrap(-5, 10) == 5
    assert index_wrap(5, 10) == 5
    assert index_wrap(15, 10) == 5


def test_uig_eval_spline_simple1():
    cs = get_cosine_spline()
    offset = 5.0
    spacing = 0.5
    uig = UniformIntGrid(np.array([-offset, -offset, -offset]), Cell(np.identity(3)*spacing), np.array([21,21,21]))

    x, y, z = np.indices(uig.shape)
    d = np.sqrt((x*spacing-offset)**2 + (y*spacing-offset)**2 + (z*spacing-offset)**2).ravel()
    data1 = cs(d)
    data2 = np.zeros(uig.shape)
    uig.eval_spline(cs, np.zeros(3), data2)
    data2 = data2.ravel()

    assert abs(data1 - data2).max() == 0.0


def test_uig_eval_spline_simple2():
    cs = get_cosine_spline()
    offset = 1.0
    spacing = 0.1
    uig = UniformIntGrid(np.array([-offset, -offset, -offset]), Cell(np.identity(3)*spacing), np.array([21,21,21]))

    x, y, z = np.indices(uig.shape)
    data1 = 0
    for ix in xrange(-1,2):
        for iy in xrange(-1,2):
            for iz in xrange(-1,2):
                d = np.sqrt((x*spacing-offset+ix*2.1)**2 + (y*spacing-offset+iy*2.1)**2 + (z*spacing-offset+iz*2.1)**2).ravel()
                data1 += cs(d)

    data2 = np.zeros(uig.shape)
    uig.eval_spline(cs, np.zeros(3), data2)
    data2 = data2.ravel()

    assert abs(data1 - data2).max() < 1e-12



def test_uig_eval_spline_with_integration():
    cs = get_cosine_spline()

    # Different integration grids
    uigs = [
        UniformIntGrid(np.array([-1.0, -1.0, -1.0]), Cell(np.identity(3)*0.1), np.array([21,21,21])),
        UniformIntGrid(np.array([-1.0, -1.0, -1.0]), Cell(np.array([[0.1, 0.1, 0.0], [0.1, 0.0, 0.1], [0.0, 0.1, 0.1]])), np.array([21,21,21])),
        UniformIntGrid(np.array([0.0, 0.0, 0.0]), Cell(np.identity(3)*0.1), np.array([21,21,21])),
        UniformIntGrid(np.array([0.0, 0.0, 0.0]), Cell(np.array([[0.1, 0.1, 0.0], [0.1, 0.0, 0.1], [0.0, 0.1, 0.1]])), np.array([21,21,21])),
        UniformIntGrid(np.array([0.0, 0.0, 0.0]), Cell(np.array([[0.1471, 0.0745, -0.010], [0.0315, -0.1403, 0.1], [0.0014, 0.147, 0.0]])), np.array([21,21,21])),
    ]

    for uig in uigs:
        # fill a 3D grid with the cosine test function
        data = np.zeros(uig.shape)
        uig.eval_spline(cs, np.zeros(3), data)

        # test the integral
        expected = 4*np.pi**2*(np.pi**2/3-2)
        assert abs(uig.integrate(data) - expected) < 1e-3


def get_simple_test_uig():
    origin = np.array([0.1, -2.0, 3.1])
    rvecs = np.array([
        [1.0, 0.1, 0.2],
        [0.1, 1.1, 0.2],
        [0.0, -0.1, 1.0],
    ])
    #origin = np.array([0.0, 0.0, 0.0])
    #rvecs = np.array([
    #    [1.0, 0.0, 0.0],
    #    [0.0, 1.0, 0.0],
    #    [0.0, 0.0, 1.0],
    #])
    shape = np.array([40, 40, 40])
    return UniformIntGrid(origin, Cell(rvecs), shape)
