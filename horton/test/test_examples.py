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


from horton.test.common import check_script, check_script_in_tmp
from horton import context


def test_example_hydrogen_ring():
    required = [context.get_fn('examples/hamiltonian/hydrogen_ring.py')]
    expected = ['ring.xyz']
    check_script_in_tmp('./hydrogen_ring.py', required, expected)


def test_example_even_tempered_li():
    check_script('./even_tempered_li.py', context.get_fn('examples/hamiltonian'))


def test_example_expectation_r():
    check_script('./expectation_r.py', context.get_fn('examples/grid'))


def test_example_rhf_water_cholesky():
    required = [context.get_fn('examples/hf_dft/rhf_water_cholesky.py'),
                context.get_fn('examples/hf_dft/water.xyz')]
    expected = ['water.h5', 'water.molden']
    check_script_in_tmp('./rhf_water_cholesky.py', required, expected)


def test_example_rhf_water_dense():
    required = [context.get_fn('examples/hf_dft/rhf_water_dense.py'),
                context.get_fn('examples/hf_dft/water.xyz')]
    expected = ['water.h5', 'water.molden']
    check_script_in_tmp('./rhf_water_dense.py', required, expected)


def test_example_rks_water_gga():
    required = [context.get_fn('examples/hf_dft/rks_water_gga.py'),
                context.get_fn('examples/hf_dft/water.xyz')]
    expected = ['water.h5']
    check_script_in_tmp('./rks_water_gga.py', required, expected)


def test_example_rks_water_hybgga():
    required = [context.get_fn('examples/hf_dft/rks_water_hybgga.py'),
                context.get_fn('examples/hf_dft/water.xyz')]
    expected = ['water.h5']
    check_script_in_tmp('./rks_water_hybgga.py', required, expected)


def test_example_rks_water_lda():
    required = [context.get_fn('examples/hf_dft/rks_water_lda.py'),
                context.get_fn('examples/hf_dft/water.xyz')]
    expected = ['water.h5']
    check_script_in_tmp('./rks_water_lda.py', required, expected)


def test_example_rks_water_numlda():
    required = [context.get_fn('examples/hf_dft/rks_water_numlda.py'),
                context.get_fn('examples/hf_dft/water.xyz')]
    expected = ['water.h5']
    check_script_in_tmp('./rks_water_numlda.py', required, expected)


def test_example_rks_water_numgga():
    required = [context.get_fn('examples/hf_dft/rks_water_numgga.py'),
                context.get_fn('examples/hf_dft/water.xyz')]
    expected = ['water.h5']
    check_script_in_tmp('./rks_water_numgga.py', required, expected)


def test_example_uhf_methyl_cholesky():
    required = [context.get_fn('examples/hf_dft/uhf_methyl_cholesky.py'),
                context.get_fn('examples/hf_dft/methyl.xyz')]
    expected = ['methyl.h5', 'methyl.molden']
    check_script_in_tmp('./uhf_methyl_cholesky.py', required, expected)


def test_example_uhf_methyl_dense():
    required = [context.get_fn('examples/hf_dft/uhf_methyl_dense.py'),
                context.get_fn('examples/hf_dft/methyl.xyz')]
    expected = ['methyl.h5', 'methyl.molden']
    check_script_in_tmp('./uhf_methyl_dense.py', required, expected)


def test_example_uks_methyl_gga():
    required = [context.get_fn('examples/hf_dft/uks_methyl_gga.py'),
                context.get_fn('examples/hf_dft/methyl.xyz')]
    expected = ['methyl.h5']
    check_script_in_tmp('./uks_methyl_gga.py', required, expected)


def test_example_uks_methyl_hybgga():
    required = [context.get_fn('examples/hf_dft/uks_methyl_hybgga.py'),
                context.get_fn('examples/hf_dft/methyl.xyz')]
    expected = ['methyl.h5']
    check_script_in_tmp('./uks_methyl_hybgga.py', required, expected)


def test_example_uks_methyl_lda():
    required = [context.get_fn('examples/hf_dft/uks_methyl_lda.py'),
                context.get_fn('examples/hf_dft/methyl.xyz')]
    expected = ['methyl.h5']
    check_script_in_tmp('./uks_methyl_lda.py', required, expected)


def test_example_uks_methyl_numlda():
    required = [context.get_fn('examples/hf_dft/uks_methyl_numlda.py'),
                context.get_fn('examples/hf_dft/methyl.xyz')]
    expected = ['methyl.h5']
    check_script_in_tmp('./uks_methyl_numlda.py', required, expected)


def test_example_uks_methyl_numgga():
    required = [context.get_fn('examples/hf_dft/uks_methyl_numgga.py'),
                context.get_fn('examples/hf_dft/methyl.xyz')]
    expected = ['methyl.h5']
    check_script_in_tmp('./uks_methyl_numgga.py', required, expected)


def test_example_wpart():
    required = [context.get_fn('examples/wpart/becke.py')]
    expected = ['charges.txt']
    check_script_in_tmp('./becke.py', required, expected)


def test_example_hf_compare():
    required = [context.get_fn('examples/hf_compare/compare.py'),
                context.get_fn('test/helium_hf_sto3g.fchk')]
    expected = []
    check_script_in_tmp('./compare.py helium_hf_sto3g.fchk', required, expected)
