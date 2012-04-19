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

from horton.cext import project_cartesian_to_pure


def test_cart_pure_s():
    work_cart = np.random.normal(0, 1, 1)
    work_pure = np.random.normal(0, 1, 1)
    project_cartesian_to_pure(work_cart, work_pure, shell_type=0, stride=1, spacing=1, count=1)
    assert abs(work_cart - work_pure).max() < 1e-10

    work_cart = np.random.normal(0, 1, 10)
    work_pure = np.random.normal(0, 1, 10)
    project_cartesian_to_pure(work_cart, work_pure, shell_type=0, stride=1, spacing=1, count=10)
    assert abs(work_cart - work_pure).max() < 1e-10


def test_cart_pure_p():
    work_cart = np.random.normal(0, 1, 3)
    work_pure = np.random.normal(0, 1, 3)
    project_cartesian_to_pure(work_cart, work_pure, shell_type=1, stride=3, spacing=1, count=1)
    assert abs(work_cart[[2,0,1]] - work_pure).max() < 1e-10

    work_cart = np.random.normal(0, 1, (10, 3))
    work_pure = np.random.normal(0, 1, (10, 3))
    project_cartesian_to_pure(work_cart.reshape(-1), work_pure.reshape(-1), shell_type=1, stride=3, spacing=1, count=10)
    assert abs(work_cart[:,[2,0,1]] - work_pure).max() < 1e-10

    work_cart = np.random.normal(0, 1, (10, 5))
    work_pure = np.random.normal(0, 1, (10, 5))
    project_cartesian_to_pure(work_cart.reshape(-1), work_pure.reshape(-1), shell_type=1, stride=5, spacing=1, count=10)
    assert abs(work_cart[:,[2,0,1]] - work_pure[:,[0,1,2]]).max() < 1e-10

    work_cart = np.random.normal(0, 1, (10, 6))
    work_pure = np.random.normal(0, 1, (10, 6))
    project_cartesian_to_pure(work_cart.reshape(-1), work_pure.reshape(-1), shell_type=1, stride=6, spacing=2, count=10)
    assert abs(work_cart[:,[4,0,2]] - work_pure[:,[0,2,4]]).max() < 1e-10

    work_cart = np.random.normal(0, 1, (3, 10))
    work_pure = np.random.normal(0, 1, (3, 10))
    project_cartesian_to_pure(work_cart.reshape(-1), work_pure.reshape(-1), shell_type=1, stride=1, spacing=10, count=10)
    assert abs(work_cart[[2,0,1],:] - work_pure).max() < 1e-10


def test_cart_pure_d():
    tf = np.array([
        [-0.5, 0, 0, -0.5, 0, 1.0],
        [0, 0, 1.0, 0, 0, 0],
        [0, 0, 0, 0, 1.0, 0],
        [0.86602540378443864676, 0, 0, -0.86602540378443864676, 0, 0],
        [0, 1.0, 0, 0, 0, 0],
    ])


    work_cart = np.random.normal(0, 1, 6)
    work_pure = np.random.normal(0, 1, 6)
    project_cartesian_to_pure(work_cart, work_pure, shell_type=2, stride=6, spacing=1, count=1)
    assert abs(np.dot(tf, work_cart) - work_pure[:5]).max() < 1e-10

    work_cart = np.random.normal(0, 1, (10, 6))
    work_pure = np.random.normal(0, 1, (10, 6))
    project_cartesian_to_pure(work_cart.reshape(-1), work_pure.reshape(-1), shell_type=2, stride=6, spacing=1, count=10)
    assert abs(np.dot(work_cart, tf.T) - work_pure[:,:5]).max() < 1e-10

    work_cart = np.random.normal(0, 1, (6, 10))
    work_pure = np.random.normal(0, 1, (6, 10))
    project_cartesian_to_pure(work_cart.reshape(-1), work_pure.reshape(-1), shell_type=2, stride=1, spacing=10, count=10)
    assert abs(np.dot(tf, work_cart) - work_pure[:5]).max() < 1e-10

def test_cart_pure_g():
    tf = np.array([
        [0.375, 0, 0, 0.21957751641341996535, 0, -0.87831006565367986142, 0, 0, 0, 0, 0.375, 0, -0.87831006565367986142, 0, 1.0],
        [0, 0, -0.89642145700079522998, 0, 0, 0, 0, -0.40089186286863657703, 0, 1.19522860933439364, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, -0.40089186286863657703, 0, 0, 0, 0, 0, 0, -0.89642145700079522998, 0, 1.19522860933439364, 0],
        [-0.5590169943749474241, 0, 0, 0, 0, 0.9819805060619657157, 0, 0, 0, 0, 0.5590169943749474241, 0, -0.9819805060619657157, 0, 0],
        [0, -0.42257712736425828875, 0, 0, 0, 0, -0.42257712736425828875, 0, 1.1338934190276816816, 0, 0, 0, 0, 0, 0],
        [0, 0, 0.790569415042094833, 0, 0, 0, 0, -1.0606601717798212866, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 1.0606601717798212866, 0, 0, 0, 0, 0, 0, -0.790569415042094833, 0, 0, 0],
        [0.73950997288745200532, 0, 0, -1.2990381056766579701, 0, 0, 0, 0, 0, 0, 0.73950997288745200532, 0, 0, 0, 0],
        [0, 1.1180339887498948482, 0, 0, 0, 0, -1.1180339887498948482, 0, 0, 0, 0, 0, 0, 0, 0],
    ])

    work_cart = np.random.normal(0, 1, 15)
    work_pure = np.random.normal(0, 1, 15)
    project_cartesian_to_pure(work_cart, work_pure, shell_type=4, stride=15, spacing=1, count=1)
    assert abs(np.dot(tf, work_cart) - work_pure[:9]).max() < 1e-10

    work_cart = np.random.normal(0, 1, (3,20))
    work_pure = np.random.normal(0, 1, (3,20))
    project_cartesian_to_pure(work_cart.reshape(-1), work_pure.reshape(-1), shell_type=4, stride=20, spacing=1, count=3)
    assert abs(np.dot(work_cart[:,:15], tf.T) - work_pure[:,:9]).max() < 1e-10
