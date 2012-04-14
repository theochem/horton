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


def test_load_operators_g09():
    lf = DenseLinalgFactory()
    eps = 1e-5
    fn = context.get_fn('test/water_sto3g_hf_g03.log')
    overlap, kinetic, nuclear_attraction, electron_repulsion = load_operators_g09(fn, lf)

    for op in overlap, kinetic, nuclear_attraction, electron_repulsion:
        assert op is not None
        assert op.nbasis == 7
        op.check_symmetry()

    assert abs(overlap.get_element(0,0) - 1.0) < eps
    assert abs(overlap.get_element(0,1) - 0.236704) < eps
    assert abs(overlap.get_element(0,2) - 0.0) < eps
    assert abs(overlap.get_element(-1,-3) - (-0.13198)) < eps

    assert abs(kinetic.get_element(2,0) - 0.0) < eps
    assert abs(kinetic.get_element(4,4) - 2.52873) < eps
    assert abs(kinetic.get_element(-1,5) - 0.00563279) < eps
    assert abs(kinetic.get_element(-1,-3) - (-0.0966161)) < eps

    assert abs(nuclear_attraction.get_element(3,3) - 9.99259) < eps
    assert abs(nuclear_attraction.get_element(-2,-1) - 1.55014) < eps
    assert abs(nuclear_attraction.get_element(2,6) - 2.76941) < eps
    assert abs(nuclear_attraction.get_element(0,3) - 0.0) < eps

    assert abs(electron_repulsion.get_element(0,0,0,0) - 4.78506575204) < eps
    assert abs(electron_repulsion.get_element(6,6,6,6) - 0.774605944194) < eps
    assert abs(electron_repulsion.get_element(6,5,0,5) - 0.0289424337101) < eps
    assert abs(electron_repulsion.get_element(5,4,0,1) - 0.00274145291476) < eps


def test_load_fchk_hf_sto3g_num():
    coordinates, numbers, basis = load_fchk(context.get_fn('test/hf_sto3g.fchk'))
    assert basis.nshell == 3
    assert basis.nbasis == 6
    assert basis._ncons.max() <= 2
    assert (basis._nexps == 3).all()
    assert len(coordinates) == len(numbers)
    assert coordinates.shape[1] == 3
    assert len(numbers) == 2


def test_load_fchk_h_sto3g_num():
    coordinates, numbers, basis = load_fchk(context.get_fn('test/h_sto3g.fchk'))
    assert basis.nshell == 1
    assert basis.nbasis == 1
    assert basis._ncons.max() <= 2
    assert (basis._nexps == 3).all()
    assert len(coordinates) == len(numbers)
    assert coordinates.shape[1] == 3
    assert len(numbers) == 1


def test_load_fchk_o2_cc_pvtz_pure_num():
    coordinates, numbers, basis = load_fchk(context.get_fn('test/o2_cc_pvtz_pure.fchk'))
    assert basis.nshell == 20
    assert basis.nbasis == 60
    assert basis._ncons.max() <= 2
    assert len(coordinates) == len(numbers)
    assert coordinates.shape[1] == 3
    assert len(numbers) == 2


def test_load_fchk_o2_cc_pvtz_cart_num():
    coordinates, numbers, basis = load_fchk(context.get_fn('test/o2_cc_pvtz_cart.fchk'))
    assert basis.nshell == 20
    assert basis.nbasis == 70
    assert basis._ncons.max() <= 2
    assert len(coordinates) == len(numbers)
    assert coordinates.shape[1] == 3
    assert len(numbers) == 2


def test_load_fchk_water_sto3g_hf():
    coordinates, numbers, basis = load_fchk(context.get_fn('test/water_sto3g_hf_g03.fchk'))
    assert basis.nshell == 4
    assert basis.nbasis == 7
    assert basis._ncons.max() == 2
    assert len(coordinates) == len(numbers)
    assert coordinates.shape[1] == 3
    assert len(numbers) == 3
