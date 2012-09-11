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


import numpy as np

from horton import *



def test_integrate_hydrogen_pair_1s():
    numbers = np.array([1, 1], int)
    centers = np.array([[0.0, 0.0, -0.5], [0.0, 0.0, 0.5]], float)
    rtf = LogRTransform(TrapezoidIntegrator1D(), 1e-3, 0.1)
    atspecs = (110, 100, rtf)

    mg = BeckeMolGrid(numbers, centers, atspecs, random_rotate=False)
    dist0 = np.sqrt(((centers[0] - mg.points)**2).sum(axis=1))
    dist1 = np.sqrt(((centers[1] - mg.points)**2).sum(axis=1))
    fn = np.exp(-2*dist0)/np.pi + np.exp(-2*dist1)/np.pi
    occupation = np.dot(fn, mg.weights)
    assert abs(occupation - 2.0) < 1e-3


def test_integrate_hydrogen_trimer_1s():
    numbers = np.array([1, 1, 1], int)
    centers = np.array([[0.0, 0.0, -0.5], [0.0, 0.0, 0.5], [0.0, 0.5, 0.0]], float)
    rtf = LogRTransform(TrapezoidIntegrator1D(), 1e-3, 0.1)
    atspecs = (110, 100, rtf)

    mg = BeckeMolGrid(numbers, centers, atspecs, random_rotate=False)
    dist0 = np.sqrt(((centers[0] - mg.points)**2).sum(axis=1))
    dist1 = np.sqrt(((centers[1] - mg.points)**2).sum(axis=1))
    dist2 = np.sqrt(((centers[2] - mg.points)**2).sum(axis=1))
    fn = np.exp(-2*dist0)/np.pi + np.exp(-2*dist1)/np.pi + np.exp(-2*dist2)/np.pi
    occupation = np.dot(fn, mg.weights)
    assert abs(occupation - 3.0) < 1e-3
