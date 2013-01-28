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


def test_cell_cubic():
    cell = Cell(np.array([[9.865, 0.0, 0.0], [0.0, 9.865, 0.0], [0.0, 0.0, 9.865]])*angstrom)
    assert cell.nvec == 3
    assert (cell.rspacings == 9.865*angstrom).all()
    assert (cell.gspacings == 1/(9.865*angstrom)).all()
    assert abs(cell.volume - (9.865*angstrom)**3) < 1e-10
    assert abs(np.dot(cell.gvecs, cell.rvecs.transpose()) - np.identity(3)).max() < 1e-5
    assert abs(np.dot(cell.gvecs.transpose(), cell.rvecs) - np.identity(3)).max() < 1e-5
    vec1 = np.array([10.0, 0.0, 5.0])*angstrom
    cell.mic(vec1)
    assert abs(vec1 - np.array([0.135, 0.0, -4.865])*angstrom).max() < 1e-10
    vec2 = np.array([10.0, 0.0, 5.0])*angstrom
    cell.add_vec(vec2, cell.to_center(vec2))
    assert abs(vec1 - vec2).max() < 1e-10
    cell.add_vec(vec1, np.array([1,2,3]))
    assert abs(vec1 - np.array([10.0, 19.73, 24.73])*angstrom).max() < 1e-10
    cell2 = Cell(-cell.rvecs)
    assert abs(cell2.volume - (9.865*angstrom)**3) < 1e-10



def test_cell_parallellogram2d():
    cell = Cell(np.array([[4.922, 0.0, 0.0], [2.462, 4.262, 0.0]])*angstrom)
    assert cell.nvec == 2
    assert abs(cell.volume - np.linalg.norm(np.cross(cell.rvecs[0], cell.rvecs[1]))) < 1e-10
    assert abs(np.dot(cell.gvecs, cell.rvecs.transpose()) - np.identity(2)).max() < 1e-5
    vec1 = np.array([10.0, 0.0, 105.0])*angstrom
    cell.mic(vec1)
    assert abs(vec1 - np.array([0.156, 0.0, 105])*angstrom).max() < 1e-3
    vec2 = np.array([10.0, 0.0, 105.0])*angstrom
    cell.add_vec(vec2, cell.to_center(vec2))
    assert abs(vec1 - vec2).max() < 1e-10
    cell.add_vec(vec1, np.array([1,2]))
    assert abs(vec1 - np.array([10.002, 8.524, 105])*angstrom).max() < 1e-3


def test_cell_1d():
    cell = Cell(np.array([[5.075, 0.187, 0.055]])*angstrom)
    assert cell.nvec == 1
    assert cell.rvecs.shape == (1, 3)
    assert cell.gvecs.shape == (1, 3)
    assert abs(cell.volume - np.linalg.norm(cell.rvecs[0])) < 1e-10
    assert abs(np.dot(cell.gvecs, cell.rvecs.transpose()) - 1) < 1e-5
    vec1 = np.array([10.0, 0.0, 105.0])*angstrom
    cell.mic(vec1)
    assert abs(vec1 - np.array([-0.15, -0.374, 104.89])*angstrom).max() < 1e-3
    vec2 = np.array([10.0, 0.0, 105.0])*angstrom
    cell.add_vec(vec2, cell.to_center(vec2))
    assert abs(vec1 - vec2).max() < 1e-10
    cell.add_vec(vec1, np.array([1]))
    assert abs(vec1 - np.array([4.925, -0.187, 104.945])*angstrom).max() < 1e-3


def test_cell_quartz():
    cell = Cell(np.array([[0.0, 0.0, 5.405222], [0.0, 4.913416, 0.0], [-4.255154, 2.456708, 0.0]])*angstrom)
    assert cell.rvecs.shape == (3, 3)
    assert cell.gvecs.shape == (3, 3)
    assert abs(cell.volume - abs(np.linalg.det(cell.rvecs))) < 1e-10
    assert abs(np.dot(cell.gvecs, cell.rvecs.transpose()) - np.identity(3)).max() < 1e-5


def test_cell_0d():
    cell = Cell()
    assert cell.nvec == 0
    assert cell.rvecs.shape == (0, 3)
    assert cell.gvecs.shape == (0, 3)
    assert cell.rspacings.shape == (0,)
    assert cell.gspacings.shape == (0,)
    vec1 = np.array([10.0, 0.0, 105.0])*angstrom
    cell.mic(vec1)
    assert abs(vec1 - np.array([10.0, 0.0, 105.0])*angstrom).max() < 1e-3
    vec2 = np.array([10.0, 0.0, 105.0])*angstrom
    cell.add_vec(vec2, cell.to_center(vec2))
    assert abs(vec1 - vec2).max() < 1e-10
    cell.add_vec(vec1, np.array([], dtype=int))
    assert abs(vec1 - np.array([10.0, 0.0, 105.0])*angstrom).max() < 1e-3
