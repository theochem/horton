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


from horton import *


def test_get_max_nbasis():
    num_contractions = np.array([1, 2, 3, 1])
    con_types = np.array([2, 1, 0, -2, 3, 0, 1])
    assert get_max_nbasis(num_contractions, con_types) == 10

    num_contractions = np.array([1, 2, 3, 1])
    con_types = np.array([1, 1, 0, -2, -2, 0, 1])
    assert get_max_nbasis(num_contractions, con_types) == 6

    num_contractions = np.array([0, 1, 1])
    con_types = np.array([1, 1])
    try:
        get_max_nbasis(num_contractions, con_types)
        assert False
    except ValueError:
        pass

    num_contractions = np.array([1, 1, 1])
    con_types = np.array([-1, 1, 2])
    try:
        get_max_nbasis(num_contractions, con_types)
        assert False
    except ValueError:
        pass

    num_contractions = np.array([1, 1, 1])
    con_types = np.array([1, 1])
    try:
        get_max_nbasis(num_contractions, con_types)
        assert False
    except ValueError:
        pass


def test_i1pow_inc_l0():
    indexes = (0,0,0)
    indexes, result = i1pow_inc(*indexes)
    assert indexes == (0,0,0)
    assert result == 0


def test_i1pow_inc_l1():
    indexes = (1,0,0)
    indexes, result = i1pow_inc(*indexes)
    assert indexes == (0,1,0)
    assert result == 1
    indexes, result = i1pow_inc(*indexes)
    assert indexes == (0,0,1)
    assert result == 1
    indexes, result = i1pow_inc(*indexes)
    assert indexes == (1,0,0)
    assert result == 0


def test_i1pow_inc_l2():
    indexes = (2,0,0)
    indexes, result = i1pow_inc(*indexes)
    assert indexes == (1,1,0)
    assert result == 1
    indexes, result = i1pow_inc(*indexes)
    assert indexes == (1,0,1)
    assert result == 1
    indexes, result = i1pow_inc(*indexes)
    assert indexes == (0,2,0)
    assert result == 1
    indexes, result = i1pow_inc(*indexes)
    assert indexes == (0,1,1)
    assert result == 1
    indexes, result = i1pow_inc(*indexes)
    assert indexes == (0,0,2)
    assert result == 1
    indexes, result = i1pow_inc(*indexes)
    assert indexes == (2,0,0)
    assert result == 0


def test_i2_pow_l0l0():
    i2p = I2Pow(0, 0, 3)
    assert i2p.fields == (0, 0, 0, 0, 0, 0, 0)
    assert i2p.inc() == False
    assert i2p.fields == (0, 0, 0, 0, 0, 0, 0)


def test_i2_pow_l1l0():
    i2p = I2Pow(1, 0, 3)
    assert i2p.fields == (1, 0, 0, 0, 0, 0, 0)
    assert i2p.inc() == True
    assert i2p.fields == (0, 1, 0, 0, 0, 0, 1)
    assert i2p.inc() == True
    assert i2p.fields == (0, 0, 1, 0, 0, 0, 2)
    assert i2p.inc() == False
    assert i2p.fields == (1, 0, 0, 0, 0, 0, 0)


def test_i2_pow_l2l1():
    i2p = I2Pow(2, 1, 10)
    assert i2p.fields == (2, 0, 0, 1, 0, 0, 0)
    assert i2p.inc() == True
    assert i2p.fields == (1, 1, 0, 1, 0, 0, 1)
    assert i2p.inc() == True
    assert i2p.fields == (1, 0, 1, 1, 0, 0, 2)
    assert i2p.inc() == True
    assert i2p.fields == (0, 2, 0, 1, 0, 0, 3)
    assert i2p.inc() == True
    assert i2p.fields == (0, 1, 1, 1, 0, 0, 4)
    assert i2p.inc() == True
    assert i2p.fields == (0, 0, 2, 1, 0, 0, 5)
    assert i2p.inc() == True
    assert i2p.fields == (2, 0, 0, 0, 1, 0, 10)
    assert i2p.inc() == True
    assert i2p.fields == (1, 1, 0, 0, 1, 0, 11)


def test_i2_pow_raise():
    try:
        i2p = I2Pow(-1,1,3)
        assert False
    except ValueError:
        pass

    try:
        i2p = I2Pow(1,-1,3)
        assert False
    except ValueError:
        pass

    try:
        i2p = I2Pow(2,1,3)
        assert False
    except ValueError:
        pass

    try:
        i2p = I2Pow(1,2,3)
        assert False
    except ValueError:
        pass
