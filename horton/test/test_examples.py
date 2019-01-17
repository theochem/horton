# -*- coding: utf-8 -*-
# HORTON: Helpful Open-source Research TOol for N-fermion systems.
# Copyright (C) 2011-2019 The HORTON Development Team
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
from nose.tools import assert_raises

from horton.test.common import check_script, check_script_in_tmp
from horton import context


def test_check_script_in_tmp():
    check_script_in_tmp('echo foo > bar', [], ['bar'])
    with assert_raises(AssertionError):
        check_script_in_tmp('echo foo > bar', [], ['egg'])
    check_script_in_tmp('echo', [context.get_fn('test/h2.xyz')], ['h2.xyz'])


@attr('rt')
def test_example_getting_started_first():
    required = [context.get_fn('examples/getting_started/first.py')]
    expected = ['water.h5']
    check_script_in_tmp('horton-regression-test.py ./first.py',
                        required, expected)


@attr('rt')
def test_example_hamiltonian_hydrogen_ring():
    required = [context.get_fn('examples/hamiltonian/hydrogen_ring.py')]
    expected = ['ring.xyz']
    check_script_in_tmp('horton-regression-test.py ./hydrogen_ring.py',
                        required, expected)


@attr('rt')
def test_example_hamiltonian_even_tempered_li():
    check_script('horton-regression-test.py ./even_tempered_li.py',
                 context.get_fn('examples/hamiltonian'))


@attr('rt')
def test_example_hamiltonian_hubbard():
    check_script('horton-regression-test.py ./hubbard.py',
                 context.get_fn('examples/hamiltonian'))


@attr('rt')
def test_example_hamiltonian_fcidump_ao():
    required = [context.get_fn('examples/hamiltonian/dump_fcidump_ao.py'),
                context.get_fn('examples/hamiltonian/load_fcidump_ao.py')]
    expected = ['hamiltonian_ao.FCIDUMP']
    check_script_in_tmp('horton-regression-test.py ./dump_fcidump_ao.py;'
                        'horton-regression-test.py ./load_fcidump_ao.py',
                        required, expected)


@attr('rt')
def test_example_hamiltonian_internal_ao():
    required = [context.get_fn('examples/hamiltonian/dump_internal_ao.py'),
                context.get_fn('examples/hamiltonian/load_internal_ao.py')]
    expected = ['hamiltonian_ao.h5']
    check_script_in_tmp('horton-regression-test.py ./dump_internal_ao.py; '
                        'horton-regression-test.py ./load_internal_ao.py',
                        required, expected)


@attr('rt')
def test_example_hamiltonian_dump_internal_ao_fcidump():
    required = [context.get_fn('examples/hamiltonian/dump_internal_ao_fcidump.py')]
    expected = ['hamiltonian_ao_fcidump.h5']
    check_script_in_tmp('horton-regression-test.py ./dump_internal_ao_fcidump.py',
                        required, expected)


@attr('rt')
def test_example_grid_expectation_r():
    check_script('horton-regression-test.py ./expectation_r.py',
                 context.get_fn('examples/grid'))


@attr('rt')
def test_example_hf_dft_rhf_h2_cholesky():
    required = [context.get_fn('examples/hf_dft/rhf_h2_cholesky.py')]
    expected = ['h2-scf.molden', 'h2-scf.h5', 'h2-hamiltonian.h5']
    check_script_in_tmp('horton-regression-test.py ./rhf_h2_cholesky.py',
                        required, expected)


@attr('rt')
def test_example_hf_dft_rhf_n2_dense():
    required = [context.get_fn('examples/hf_dft/rhf_n2_dense.py')]
    expected = ['n2-scf.molden', 'n2-scf.h5', 'n2-cas8-8.FCIDUMP',
                'n2-hamiltonian-cas8-8.h5', 'n2.FCIDUMP', 'n2-hamiltonian.h5']
    check_script_in_tmp('horton-regression-test.py ./rhf_n2_dense.py',
                        required, expected)


@attr('rt')
def test_example_hf_dft_rhf_water_cholesky():
    required = [context.get_fn('examples/hf_dft/rhf_water_cholesky.py')]
    expected = ['water.h5', 'water.molden']
    check_script_in_tmp('horton-regression-test.py ./rhf_water_cholesky.py',
                        required, expected)


@attr('rt')
def test_example_hf_dft_rhf_water_dense():
    required = [context.get_fn('examples/hf_dft/rhf_water_dense.py')]
    expected = ['water-scf.h5', 'water-scf.molden', 'water.FCIDUMP',
                'water-hamiltonian.h5']
    check_script_in_tmp('horton-regression-test.py ./rhf_water_dense.py',
                        required, expected)


@attr('rt')
def test_example_hf_dft_rks_water_gga():
    required = [context.get_fn('examples/hf_dft/rks_water_gga.py')]
    expected = ['water.h5']
    check_script_in_tmp('horton-regression-test.py ./rks_water_gga.py',
                        required, expected)


@attr('rt')
def test_example_hf_dft_rks_water_hybgga():
    required = [context.get_fn('examples/hf_dft/rks_water_hybgga.py')]
    expected = ['water.h5']
    check_script_in_tmp('horton-regression-test.py ./rks_water_hybgga.py',
                        required, expected)


@attr('rt')
def test_example_hf_dft_rks_water_hybmgga():
    required = [context.get_fn('examples/hf_dft/rks_water_hybmgga.py')]
    expected = ['water.h5']
    check_script_in_tmp('horton-regression-test.py ./rks_water_hybmgga.py',
                        required, expected)


@attr('rt')
def test_example_hf_dft_rks_water_lda():
    required = [context.get_fn('examples/hf_dft/rks_water_lda.py')]
    expected = ['water.h5']
    check_script_in_tmp('horton-regression-test.py ./rks_water_lda.py',
                        required, expected)


@attr('rt')
def test_example_hf_dft_rks_water_mgga():
    required = [context.get_fn('examples/hf_dft/rks_water_mgga.py')]
    expected = ['water.h5']
    check_script_in_tmp('horton-regression-test.py ./rks_water_mgga.py',
                        required, expected)


@attr('rt')
def test_example_hf_dft_rks_water_numgga():
    required = [context.get_fn('examples/hf_dft/rks_water_numgga.py')]
    expected = ['water.h5']
    check_script_in_tmp('horton-regression-test.py ./rks_water_numgga.py',
                        required, expected)


@attr('rt')
def test_example_hdf_dft_rks_water_numlda():
    required = [context.get_fn('examples/hf_dft/rks_water_numlda.py')]
    expected = ['water.h5']
    check_script_in_tmp('horton-regression-test.py ./rks_water_numlda.py',
                        required, expected)


@attr('rt')
def test_example_hf_dft_uhf_methyl_cholesky():
    required = [context.get_fn('examples/hf_dft/uhf_methyl_cholesky.py')]
    expected = ['methyl.h5', 'methyl.molden']
    check_script_in_tmp('horton-regression-test.py ./uhf_methyl_cholesky.py',
                        required, expected)


@attr('rt')
def test_example_hf_dft_uhf_methyl_dense():
    required = [context.get_fn('examples/hf_dft/uhf_methyl_dense.py')]
    expected = ['methyl.h5', 'methyl.molden']
    check_script_in_tmp('horton-regression-test.py ./uhf_methyl_dense.py',
                        required, expected)


@attr('rt')
def test_example_hf_dft_uks_methyl_gga():
    required = [context.get_fn('examples/hf_dft/uks_methyl_gga.py')]
    expected = ['methyl.h5']
    check_script_in_tmp('horton-regression-test.py ./uks_methyl_gga.py',
                        required, expected)


@attr('rt')
def test_example_hf_dft_uks_methyl_hybgga():
    required = [context.get_fn('examples/hf_dft/uks_methyl_hybgga.py')]
    expected = ['methyl.h5']
    check_script_in_tmp('horton-regression-test.py ./uks_methyl_hybgga.py',
                        required, expected)


@attr('rt')
def test_example_hf_dft_uks_methyl_hybmgga():
    required = [context.get_fn('examples/hf_dft/uks_methyl_hybmgga.py')]
    expected = ['methyl.h5']
    check_script_in_tmp('horton-regression-test.py ./uks_methyl_hybmgga.py',
                        required, expected)


@attr('rt')
def test_example_hf_dft_uks_methyl_lda():
    required = [context.get_fn('examples/hf_dft/uks_methyl_lda.py')]
    expected = ['methyl.h5']
    check_script_in_tmp('horton-regression-test.py ./uks_methyl_lda.py',
                        required, expected)


@attr('rt')
def test_example_hf_dft_uks_methyl_mgga():
    required = [context.get_fn('examples/hf_dft/uks_methyl_mgga.py')]
    expected = ['methyl.h5']
    check_script_in_tmp('horton-regression-test.py ./uks_methyl_mgga.py',
                        required, expected)


@attr('rt')
def test_example_hf_dft_uks_methyl_numgga():
    required = [context.get_fn('examples/hf_dft/uks_methyl_numgga.py')]
    expected = ['methyl.h5']
    check_script_in_tmp('horton-regression-test.py ./uks_methyl_numgga.py',
                        required, expected)


@attr('rt')
def test_example_hf_dft_uks_methyl_numlda():
    required = [context.get_fn('examples/hf_dft/uks_methyl_numlda.py')]
    expected = ['methyl.h5']
    check_script_in_tmp('horton-regression-test.py ./uks_methyl_numlda.py',
                        required, expected)


@attr('rt')
def test_example_wpart_becke():
    required = [context.get_fn('examples/wpart/becke.py')]
    expected = ['charges.txt']
    check_script_in_tmp('horton-regression-test.py ./becke.py',
                        required, expected)


@attr('rt')
def test_example_hf_compare():
    required = [context.get_fn('examples/hf_compare/compare.py'),
                context.get_fn('examples/hf_compare/compare_rt.py')]
    expected = ['compare.log']
    cmd = './compare.py %s > compare.log' % context.get_fn('test/helium_hf_sto3g.fchk')
    cmd += '; horton-regression-test.py ./compare_rt.py'
    check_script_in_tmp(cmd, required, expected)
