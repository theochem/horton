# -*- coding: utf-8 -*-
# HORTON: Helpful Open-source Research TOol for N-fermion systems.
# Copyright (C) 2011-2015 The HORTON Development Team
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
#--
#pylint: skip-file


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


@attr('slow')
def test_example_ap1rog_hubbard():
    required = [context.get_fn('examples/ap1rog/hubbard.py')]
    expected = ['checkpoint.h5']
    check_script_in_tmp('./hubbard.py', required, expected)


@attr('slow')
def test_example_ap1rog_extham_n2():
    required = [context.get_fn('examples/ap1rog/external_hamiltonian_n2_dense.py'),
                context.get_fn('examples/hf_dft/rhf_n2_dense.py')]
    expected = ['checkpoint.h5', 'n2-scf.molden', 'n2-scf.h5',
                'n2-cas8-8.FCIDUMP', 'n2-hamiltonian-cas8-8.h5', 'n2.FCIDUMP',
                'n2-hamiltonian.h5']
    check_script_in_tmp('./rhf_n2_dense.py; ./external_hamiltonian_n2_dense.py', required, expected)


@attr('slow')
def test_example_ap1rog_extham_h2():
    required = [context.get_fn('examples/ap1rog/external_hamiltonian_h2_cholesky.py'),
                context.get_fn('examples/hf_dft/rhf_h2_cholesky.py')]
    expected = ['checkpoint.h5', 'h2-scf.molden', 'h2-scf.h5',
                'h2-hamiltonian.h5']
    check_script_in_tmp('./rhf_h2_cholesky.py; ./external_hamiltonian_h2_cholesky.py', required, expected)


@attr('slow')
def test_example_ap1rog_extham_water():
    required = [context.get_fn('examples/ap1rog/external_hamiltonian_water_dense.py'),
                context.get_fn('examples/hf_dft/rhf_water_dense.py')]
    expected = ['checkpoint.h5', 'water-scf.molden', 'water-scf.h5',
                'water-hamiltonian.h5']
    check_script_in_tmp('./rhf_water_dense.py; ./external_hamiltonian_water_dense.py', required, expected)


@attr('slow')
def test_example_ap1rog_h2_cholesky_321g():
    required = [context.get_fn('examples/ap1rog/h2_cholesky_3-21g.py')]
    expected = ['checkpoint.h5']
    check_script_in_tmp('./h2_cholesky_3-21g.py', required, expected)


@attr('slow')
def test_example_ap1rog_h2_cholesky_aug_cc_pvdz():
    required = [context.get_fn('examples/ap1rog/h2_cholesky_aug-cc-pvdz.py')]
    expected = ['checkpoint.h5']
    check_script_in_tmp('./h2_cholesky_aug-cc-pvdz.py', required, expected)


@attr('slow')
def test_example_ap1rog_h2_cholesky_cc_pvdz():
    required = [context.get_fn('examples/ap1rog/h2_cholesky_cc-pvdz.py')]
    expected = ['checkpoint.h5']
    check_script_in_tmp('./h2_cholesky_cc-pvdz.py', required, expected)


@attr('slow')
def test_example_ap1rog_h2_dense_321g():
    required = [context.get_fn('examples/ap1rog/h2_dense_3-21g.py')]
    expected = ['checkpoint.h5']
    check_script_in_tmp('./h2_dense_3-21g.py', required, expected)


@attr('slow')
def test_example_ap1rog_h2_dense_aug_cc_pvdz():
    required = [context.get_fn('examples/ap1rog/h2_dense_aug-cc-pvdz.py')]
    expected = ['checkpoint.h5']
    check_script_in_tmp('./h2_dense_aug-cc-pvdz.py', required, expected)


@attr('slow')
def test_example_ap1rog_water_cholesky_321g():
    required = [context.get_fn('examples/ap1rog/water_cholesky_3-21g.py')]
    expected = ['checkpoint.h5']
    check_script_in_tmp('./water_cholesky_3-21g.py', required, expected)


@attr('slow')
def test_example_ap1rog_water_cholesky_321g_restart():
    required = [context.get_fn('examples/ap1rog/water_cholesky_3-21g.py'),
                context.get_fn('examples/ap1rog/restart_water_cholesky_3-21g.py')]
    expected = ['checkpoint.h5']
    check_script_in_tmp('./water_cholesky_3-21g.py; ./restart_water_cholesky_3-21g.py', required, expected)


@attr('slow')
def test_example_ap1rog_water_cholesky_sto3g():
    required = [context.get_fn('examples/ap1rog/water_cholesky_sto-3g.py')]
    expected = ['checkpoint.h5']
    check_script_in_tmp('./water_cholesky_sto-3g.py', required, expected)


@attr('slow')
def test_example_ap1rog_water_cholesky_cc_pvdz():
    required = [context.get_fn('examples/ap1rog/water_cholesky_cc-pvdz.py')]
    expected = ['checkpoint.h5']
    check_script_in_tmp('./water_cholesky_cc-pvdz.py', required, expected)


@attr('slow')
def test_example_ap1rog_water_dense_321g():
    required = [context.get_fn('examples/ap1rog/water_dense_3-21g.py')]
    expected = ['checkpoint.h5']
    check_script_in_tmp('./water_dense_3-21g.py', required, expected)


@attr('slow')
def test_example_ap1rog_water_dense_sto3g():
    required = [context.get_fn('examples/ap1rog/water_dense_sto-3g.py')]
    expected = ['checkpoint.h5']
    check_script_in_tmp('./water_dense_sto-3g.py', required, expected)


@attr('slow')
def test_example_ap1rog_water_dense_cc_pvdz():
    required = [context.get_fn('examples/ap1rog/water_dense_cc-pvdz.py')]
    expected = ['checkpoint.h5']
    check_script_in_tmp('./water_dense_cc-pvdz.py', required, expected)


@attr('slow')
def test_example_ap1rog_water_default():
    required = [context.get_fn('examples/ap1rog/water_default.py')]
    expected = ['checkpoint.h5']
    check_script_in_tmp('./water_default.py', required, expected)


@attr('slow')
def test_example_ap1rog_water_minimal():
    required = [context.get_fn('examples/ap1rog/water_minimal.py')]
    expected = ['checkpoint.h5']
    check_script_in_tmp('./water_minimal.py', required, expected)


def test_example_localization():
    check_script('./water_pm.py', context.get_fn('examples/localization'))


def test_example_pt_mp2_h2_431g():
    check_script('./mp2_h2_4-31g.py', context.get_fn('examples/perturbation_theory'))


def test_example_pt_mp2_h2_cc_pvdz():
    check_script('./mp2_h2_cc-pvdz.py', context.get_fn('examples/perturbation_theory'))


def test_example_pt_mp2_h2_def2_svpd():
    check_script('./mp2_h2_def2-svpd.py', context.get_fn('examples/perturbation_theory'))


def test_example_pt_mp2_water_431g():
    check_script('./mp2_water_4-31g.py', context.get_fn('examples/perturbation_theory'))


def test_example_pt_mp2_water_cc_pvdz():
    check_script('./mp2_water_cc-pvdz.py', context.get_fn('examples/perturbation_theory'))


@attr('slow')
def test_example_pt_mp2_water_def2_svpd():
    check_script('./mp2_water_def2-svpd.py', context.get_fn('examples/perturbation_theory'))


@attr('slow')
def test_example_pta_h2_431g():
    required = [context.get_fn('examples/perturbation_theory/pta_h2_4-31g.py')]
    expected = ['checkpoint.h5']
    check_script_in_tmp('./pta_h2_4-31g.py', required, expected)


@attr('slow')
def test_example_pta_h2_cc_pvdz():
    required = [context.get_fn('examples/perturbation_theory/pta_h2_cc-pvdz.py')]
    expected = ['checkpoint.h5']
    check_script_in_tmp('./pta_h2_cc-pvdz.py', required, expected)


@attr('slow')
def test_example_pta_h2_def2_svpd():
    required = [context.get_fn('examples/perturbation_theory/pta_h2_def2-svpd.py')]
    expected = ['checkpoint.h5']
    check_script_in_tmp('./pta_h2_def2-svpd.py', required, expected)


@attr('slow')
def test_example_pta_water_431g():
    required = [context.get_fn('examples/perturbation_theory/pta_water_4-31g.py')]
    expected = ['checkpoint.h5']
    check_script_in_tmp('./pta_water_4-31g.py', required, expected)


@attr('slow')
def test_example_pta_water_cc_pvdz():
    required = [context.get_fn('examples/perturbation_theory/pta_water_cc-pvdz.py')]
    expected = ['checkpoint.h5']
    check_script_in_tmp('./pta_water_cc-pvdz.py', required, expected)


@attr('slow')
def test_example_pta_water_def2_svpd():
    required = [context.get_fn('examples/perturbation_theory/pta_water_def2-svpd.py')]
    expected = ['checkpoint.h5']
    check_script_in_tmp('./pta_water_def2-svpd.py', required, expected)


@attr('slow')
def test_example_ptb_h2_431g():
    required = [context.get_fn('examples/perturbation_theory/ptb_h2_4-31g.py')]
    expected = ['checkpoint.h5']
    check_script_in_tmp('./ptb_h2_4-31g.py', required, expected)


@attr('slow')
def test_example_ptb_h2_cc_pvdz():
    required = [context.get_fn('examples/perturbation_theory/ptb_h2_cc-pvdz.py')]
    expected = ['checkpoint.h5']
    check_script_in_tmp('./ptb_h2_cc-pvdz.py', required, expected)


@attr('slow')
def test_example_ptb_h2_def2_svpd():
    required = [context.get_fn('examples/perturbation_theory/ptb_h2_def2-svpd.py')]
    expected = ['checkpoint.h5']
    check_script_in_tmp('./ptb_h2_def2-svpd.py', required, expected)


@attr('slow')
def test_example_ptb_water_431g():
    required = [context.get_fn('examples/perturbation_theory/ptb_water_4-31g.py')]
    expected = ['checkpoint.h5']
    check_script_in_tmp('./ptb_water_4-31g.py', required, expected)


@attr('slow')
def test_example_ptb_water_cc_pvdz():
    required = [context.get_fn('examples/perturbation_theory/ptb_water_cc-pvdz.py')]
    expected = ['checkpoint.h5']
    check_script_in_tmp('./ptb_water_cc-pvdz.py', required, expected)


@attr('slow')
def test_example_ptb_water_def2_svpd():
    required = [context.get_fn('examples/perturbation_theory/ptb_water_def2-svpd.py')]
    expected = ['checkpoint.h5']
    check_script_in_tmp('./ptb_water_def2-svpd.py', required, expected)


@attr('slow')
def test_example_oe_water():
    required = [context.get_fn('examples/orbital_entanglement/water.py')]
    expected = ['i12.dat', 'checkpoint.h5', 's1.dat', 'orbital_entanglement.png']
    check_script_in_tmp('./water.py; horton-entanglement.py 0.001; '
                        'horton-entanglement.py 0.001 1 3;'
                        'horton-entanglement.py 0.0001 10', required, expected)


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


def test_example_hf_compare():
    required = [context.get_fn('examples/hf_compare/compare.py')]
    expected = []
    check_script_in_tmp('./compare.py %s' % context.get_fn('test/helium_hf_sto3g.fchk'), required, expected)
