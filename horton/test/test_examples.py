# -*- coding: utf-8 -*-
# Horton is a development platform for electronic structure methods.
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
#pylint: skip-file


from horton.test.common import check_script
from horton import context


def test_example_hydrogen_ring():
    check_script('./hydrogen_ring.py', context.get_fn('examples/hamiltonian'))


def test_example_even_tempered_li():
    check_script('./even_tempered_li.py', context.get_fn('examples/hamiltonian'))


def test_example_rhf_cr2_cholesky():
    check_script('./rhf_cr2_cholesky.py', context.get_fn('examples/hf_dft'))


def test_example_rhf_water_cholesky():
    check_script('./rhf_water_cholesky.py', context.get_fn('examples/hf_dft'))


def test_example_rhf_water_dense():
    check_script('./rhf_water_dense.py', context.get_fn('examples/hf_dft'))


def test_example_rks_water_gga():
    check_script('./rks_water_gga.py', context.get_fn('examples/hf_dft'))


def test_example_rks_water_hybgga():
    check_script('./rks_water_hybgga.py', context.get_fn('examples/hf_dft'))


def test_example_rks_water_lda():
    check_script('./rks_water_lda.py', context.get_fn('examples/hf_dft'))


def test_example_rks_water_numlda():
    check_script('./rks_water_numlda.py', context.get_fn('examples/hf_dft'))


def test_example_uhf_methyl_cholesky():
    check_script('./uhf_methyl_cholesky.py', context.get_fn('examples/hf_dft'))


def test_example_uhf_methyl_dense():
    check_script('./uhf_methyl_dense.py', context.get_fn('examples/hf_dft'))


def test_example_uks_methyl_gga():
    check_script('./uks_methyl_gga.py', context.get_fn('examples/hf_dft'))


def test_example_uks_methyl_hybgga():
    check_script('./uks_methyl_hybgga.py', context.get_fn('examples/hf_dft'))


def test_example_uks_methyl_lda():
    check_script('./uks_methyl_lda.py', context.get_fn('examples/hf_dft'))


def test_example_uks_methyl_numlda():
    check_script('./uks_methyl_numlda.py', context.get_fn('examples/hf_dft'))


def test_example_wpart():
    check_script('./becke.py', context.get_fn('examples/wpart'))


def test_example_hf_compare():
    fn_fchk = context.get_fn('test/helium_hf_sto3g.fchk')
    check_script('./compare.py %s' % fn_fchk, context.get_fn('examples/hf_compare'))
