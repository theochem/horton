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


from nose.plugins.attrib import attr

from horton.test.common import check_script, check_script_in_tmp
from horton import context


def test_example_first():
    required = [context.get_fn('examples/getting_started/first.py')]
    expected = ['water.h5']
    check_script_in_tmp('./first.py', required, expected)


def test_example_hydrogen_ring():
    required = [context.get_fn('examples/hamiltonian/hydrogen_ring.py')]
    expected = ['ring.xyz']
    check_script_in_tmp('./hydrogen_ring.py', required, expected)


def test_example_even_tempered_li():
    check_script('./even_tempered_li.py', context.get_fn('examples/hamiltonian'))


def test_example_hubbard():
    check_script('./hubbard.py', context.get_fn('examples/hamiltonian'))


@attr('slow')
def test_hamiltonian_fcidump_ao():
    required = [context.get_fn('examples/hamiltonian/dump_fcidump_ao.py'),
                context.get_fn('examples/hamiltonian/load_fcidump_ao.py')]
    expected = ['hamiltonian_ao.FCIDUMP']
    check_script_in_tmp('./dump_fcidump_ao.py; ./load_fcidump_ao.py', required, expected)


def test_hamiltonian_internal_ao():
    required = [context.get_fn('examples/hamiltonian/dump_internal_ao.py'),
                context.get_fn('examples/hamiltonian/load_internal_ao.py')]
    expected = ['hamiltonian_ao.h5']
    check_script_in_tmp('./dump_internal_ao.py; ./load_internal_ao.py', required, expected)


def test_hamiltonian_internal_ao_fcidump():
    required = [context.get_fn('examples/hamiltonian/dump_internal_ao_fcidump.py')]
    expected = ['hamiltonian_ao_fcidump.h5']
    check_script_in_tmp('./dump_internal_ao_fcidump.py', required, expected)


def test_example_expectation_r():
    check_script('./expectation_r.py', context.get_fn('examples/grid'))


def test_example_rhf_water_cholesky():
    required = [context.get_fn('examples/hf_dft/rhf_water_cholesky.py')]
    expected = ['water.h5', 'water.molden']
    check_script_in_tmp('./rhf_water_cholesky.py', required, expected)


@attr('slow')
def test_example_rhf_n2_dense():
    required = [context.get_fn('examples/hf_dft/rhf_n2_dense.py')]
    expected = ['n2-scf.molden', 'n2-scf.h5', 'n2-cas8-8.FCIDUMP',
                'n2-hamiltonian-cas8-8.h5', 'n2.FCIDUMP', 'n2-hamiltonian.h5']
    check_script_in_tmp('./rhf_n2_dense.py', required, expected)


def test_example_rhf_h2_cholesky():
    required = [context.get_fn('examples/hf_dft/rhf_h2_cholesky.py')]
    expected = ['h2-scf.molden', 'h2-scf.h5', 'h2-hamiltonian.h5']
    check_script_in_tmp('./rhf_h2_cholesky.py', required, expected)


def test_example_rhf_water_dense():
    required = [context.get_fn('examples/hf_dft/rhf_water_dense.py')]
    expected = ['water-scf.h5', 'water-scf.molden', 'water.FCIDUMP',
                'water-hamiltonian.h5']
    check_script_in_tmp('./rhf_water_dense.py', required, expected)


@attr('slow')
def test_example_rks_water_gga():
    required = [context.get_fn('examples/hf_dft/rks_water_gga.py')]
    expected = ['water.h5']
    check_script_in_tmp('./rks_water_gga.py', required, expected)


@attr('slow')
def test_example_rks_water_hybgga():
    required = [context.get_fn('examples/hf_dft/rks_water_hybgga.py')]
    expected = ['water.h5']
    check_script_in_tmp('./rks_water_hybgga.py', required, expected)


@attr('slow')
def test_example_rks_water_lda():
    required = [context.get_fn('examples/hf_dft/rks_water_lda.py')]
    expected = ['water.h5']
    check_script_in_tmp('./rks_water_lda.py', required, expected)


@attr('slow')
def test_example_rks_water_numlda():
    required = [context.get_fn('examples/hf_dft/rks_water_numlda.py')]
    expected = ['water.h5']
    check_script_in_tmp('./rks_water_numlda.py', required, expected)


@attr('slow')
def test_example_rks_water_numgga():
    required = [context.get_fn('examples/hf_dft/rks_water_numgga.py')]
    expected = ['water.h5']
    check_script_in_tmp('./rks_water_numgga.py', required, expected)


def test_example_uhf_methyl_cholesky():
    required = [context.get_fn('examples/hf_dft/uhf_methyl_cholesky.py')]
    expected = ['methyl.h5', 'methyl.molden']
    check_script_in_tmp('./uhf_methyl_cholesky.py', required, expected)


def test_example_uhf_methyl_dense():
    required = [context.get_fn('examples/hf_dft/uhf_methyl_dense.py')]
    expected = ['methyl.h5', 'methyl.molden']
    check_script_in_tmp('./uhf_methyl_dense.py', required, expected)


@attr('slow')
def test_example_uks_methyl_gga():
    required = [context.get_fn('examples/hf_dft/uks_methyl_gga.py')]
    expected = ['methyl.h5']
    check_script_in_tmp('./uks_methyl_gga.py', required, expected)


@attr('slow')
def test_example_uks_methyl_hybgga():
    required = [context.get_fn('examples/hf_dft/uks_methyl_hybgga.py')]
    expected = ['methyl.h5']
    check_script_in_tmp('./uks_methyl_hybgga.py', required, expected)


@attr('slow')
def test_example_uks_methyl_lda():
    required = [context.get_fn('examples/hf_dft/uks_methyl_lda.py')]
    expected = ['methyl.h5']
    check_script_in_tmp('./uks_methyl_lda.py', required, expected)


@attr('slow')
def test_example_uks_methyl_numlda():
    required = [context.get_fn('examples/hf_dft/uks_methyl_numlda.py')]
    expected = ['methyl.h5']
    check_script_in_tmp('./uks_methyl_numlda.py', required, expected)


@attr('slow')
def test_example_uks_methyl_numgga():
    required = [context.get_fn('examples/hf_dft/uks_methyl_numgga.py')]
    expected = ['methyl.h5']
    check_script_in_tmp('./uks_methyl_numgga.py', required, expected)


def test_example_wpart():
    required = [context.get_fn('examples/wpart/becke.py')]
    expected = ['charges.txt']
    check_script_in_tmp('./becke.py', required, expected)


@attr('slow')
def test_example_hf_compare():
    required = [context.get_fn('examples/hf_compare/compare.py')]
    expected = []
    check_script_in_tmp('./compare.py %s' % context.get_fn('test/helium_hf_sto3g.fchk'), required, expected)
