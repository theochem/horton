# -*- coding: utf-8 -*-
# Horton is a Density Functional Theory program.
# Copyright (C) 2011 Toon Verstraelen <Toon.Verstraelen@UGent.be>, ...
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


from horton import context, GaussianOrbBasis, get_con_nbasis

from molmod.io import FCHKFile


def test_hf_sto3g_num():
    fchk = FCHKFile(context.get_fn('test/hf_sto3g.fchk'))
    basis = GaussianOrbBasis.from_fchk(fchk)
    assert len(basis.centers) == 2
    assert basis.nshell == 3
    assert basis.nbasis == 6
    assert basis._num_contractions.max() <= 2
    assert (basis._num_exponents == 3).all()


def test_h_sto3g_num():
    fchk = FCHKFile(context.get_fn('test/h_sto3g.fchk'))
    basis = GaussianOrbBasis.from_fchk(fchk)
    assert len(basis.centers) == 1
    assert basis.nshell == 1
    assert basis.nbasis == 1
    assert basis._num_contractions.max() <= 2
    assert (basis._num_exponents == 3).all()


def test_o2_cc_pvtz_pure_num():
    fchk = FCHKFile(context.get_fn('test/o2_cc_pvtz_pure.fchk'))
    basis = GaussianOrbBasis.from_fchk(fchk)
    assert len(basis.centers)==2
    assert basis.nshell == 20
    assert basis.nbasis == 60
    assert basis._num_contractions.max() <= 2


def test_o2_cc_pvtz_cart_num():
    fchk = FCHKFile(context.get_fn('test/o2_cc_pvtz_cart.fchk'))
    basis = GaussianOrbBasis.from_fchk(fchk)
    assert len(basis.centers)==2
    assert basis.nshell == 20
    assert basis.nbasis == 70
    assert basis._num_contractions.max() <= 2


def test_con_nbasis():
    assert get_con_nbasis(-3) == 7
    assert get_con_nbasis(-2) == 5
    assert get_con_nbasis( 0) == 1
    assert get_con_nbasis( 1) == 3
    assert get_con_nbasis( 2) == 6
    assert get_con_nbasis( 3) == 10
    try:
        get_con_nbasis(-1)
        assert False
    except ValueError:
        pass
