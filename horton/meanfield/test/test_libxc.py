# -*- coding: utf-8 -*-
# HORTON: Helpful Open-source Research TOol for N-fermion systems.
# Copyright (C) 2011-2022 The HORTON Development Team
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


from nose.tools import assert_raises
from nose.plugins.skip import SkipTest

from horton import *  # pylint: disable=wildcard-import,unused-wildcard-import
from horton.meanfield.test.common import check_interpolation, \
    check_dot_hessian, check_dot_hessian_polynomial, check_dot_hessian_cache


def setup_gga_cs(name):
    """Prepare datastructures for R-GGA calculation in CO."""
    fn_fchk = context.get_fn('test/co_pbe_sto3g.fchk')
    mol = IOData.from_file(fn_fchk)
    mol.dm_alpha = mol.orb_alpha.to_dm()

    grid = BeckeMolGrid(mol.coordinates, mol.numbers, mol.pseudo_numbers, random_rotate=False)
    olp = mol.obasis.compute_overlap()
    kin = mol.obasis.compute_kinetic()
    na = mol.obasis.compute_nuclear_attraction(mol.coordinates, mol.pseudo_numbers)
    terms = [
        RGridGroup(mol.obasis, grid, [
            RLibXCGGA(name),
        ]),
    ]
    ham = REffHam(terms)
    return mol, olp, kin, na, ham


def test_cubic_interpolation_c_pbe_cs():
    mol, olp, kin, na, ham = setup_gga_cs('c_pbe')
    check_interpolation(ham, olp, kin, na, [mol.orb_alpha])


def test_dot_hessian_c_pbe_cs():
    mol, _olp, _kin, _na, ham = setup_gga_cs('c_pbe')
    check_dot_hessian(ham, mol.dm_alpha)


def test_dot_hessian_c_pbe_cs_polynomial():
    mol, olp, kin, na, ham = setup_gga_cs('c_pbe')
    check_dot_hessian_polynomial(olp, kin+na, ham, [mol.orb_alpha], is_hf=False, extent=0.00005)


def test_dot_hessian_c_pbe_cs_cache():
    mol, _olp, _kin, _na, ham = setup_gga_cs('c_pbe')
    check_dot_hessian_cache(ham, mol.dm_alpha)


def test_cubic_interpolation_x_pbe_cs():
    mol, olp, kin, na, ham = setup_gga_cs('x_pbe')
    check_interpolation(ham, olp, kin, na, [mol.orb_alpha])


def test_dot_hessian_x_pbe_cs():
    mol, _olp, _kin, _na, ham = setup_gga_cs('x_pbe')
    check_dot_hessian(ham, mol.dm_alpha)


def test_dot_hessian_x_pbe_cs_polynomial():
    mol, olp, kin, na, ham = setup_gga_cs('x_pbe')
    check_dot_hessian_polynomial(olp, kin+na, ham, [mol.orb_alpha], is_hf=False, extent=0.00005)


def test_dot_hessian_x_pbe_cs_cache():
    mol, _olp, _kin, _na, ham = setup_gga_cs('x_pbe')
    check_dot_hessian_cache(ham, mol.dm_alpha)


def setup_hfs_cs():
    """Prepare datastructures for R-HFS (x-functional-only) calculation on CO."""
    fn_fchk = context.get_fn('test/co_pbe_sto3g.fchk')
    mol = IOData.from_file(fn_fchk)
    mol.dm_alpha = mol.orb_alpha.to_dm()

    grid = BeckeMolGrid(mol.coordinates, mol.numbers, mol.pseudo_numbers,
                        'coarse', random_rotate=False)
    olp = mol.obasis.compute_overlap()
    kin = mol.obasis.compute_kinetic()
    na = mol.obasis.compute_nuclear_attraction(mol.coordinates, mol.pseudo_numbers)
    terms = [
        RGridGroup(mol.obasis, grid, [
            RLibXCLDA('x'),
        ]),
    ]
    ham = REffHam(terms)

    return mol, olp, kin, na, ham


def test_cubic_interpolation_hfs_cs():
    mol, olp, kin, na, ham = setup_hfs_cs()
    check_interpolation(ham, olp, kin, na, [mol.orb_alpha])


def test_dot_hessian_hfs_cs():
    mol, _olp, _kin, _na, ham = setup_hfs_cs()
    check_dot_hessian(ham, mol.dm_alpha)


def test_dot_hessian_hfs_cs_polynomial():
    mol, olp, kin, na, ham = setup_hfs_cs()
    check_dot_hessian_polynomial(olp, kin+na, ham, [mol.orb_alpha], is_hf=False, extent=0.00005)


def test_dot_hessian_hfs_cs_cache():
    mol, _olp, _kin, _na, ham = setup_hfs_cs()
    check_dot_hessian_cache(ham, mol.dm_alpha)


def setup_x_pbe_c_vwn_cs():
    """Setup datastructure for mixed GGA+LDA calculation."""
    # mixing of GGA and LDA
    fn_fchk = context.get_fn('test/water_hfs_321g.fchk')
    mol = IOData.from_file(fn_fchk)
    mol.dm_alpha = mol.orb_alpha.to_dm()

    grid = BeckeMolGrid(mol.coordinates, mol.numbers, mol.pseudo_numbers, random_rotate=False)
    olp = mol.obasis.compute_overlap()
    kin = mol.obasis.compute_kinetic()
    na = mol.obasis.compute_nuclear_attraction(mol.coordinates, mol.pseudo_numbers)
    terms = [
        RGridGroup(mol.obasis, grid, [
            RLibXCGGA('x_pbe'),
            RLibXCLDA('c_vwn'),
        ]),
    ]
    ham = REffHam(terms)
    return mol, olp, kin, na, ham


def test_cubic_interpolation_x_pbe_c_vwn_cs():
    mol, olp, kin, na, ham = setup_x_pbe_c_vwn_cs()
    check_interpolation(ham, olp, kin, na, [mol.orb_alpha])


def test_dot_hessian_x_pbe_c_vwn_cs():
    mol, _olp, _kin, _na, ham = setup_x_pbe_c_vwn_cs()
    check_dot_hessian(ham, mol.dm_alpha)


def test_dot_hessian_x_pbe_c_vwn_cs_polynomial():
    mol, olp, kin, na, ham = setup_x_pbe_c_vwn_cs()
    check_dot_hessian_polynomial(olp, kin+na, ham, [mol.orb_alpha], is_hf=False, extent=0.00005)


def test_dot_hessian_x_pbe_c_vwn_cs_cache():
    mol, _olp, _kin, _na, ham = setup_x_pbe_c_vwn_cs()
    check_dot_hessian_cache(ham, mol.dm_alpha)


def setup_c_vwn_cs():
    """Prepare datastructures for R-VWN (c-functional-only) calculation on water."""
    fn_fchk = context.get_fn('test/water_hfs_321g.fchk')
    mol = IOData.from_file(fn_fchk)
    mol.dm_alpha = mol.orb_alpha.to_dm()

    grid = BeckeMolGrid(mol.coordinates, mol.numbers, mol.pseudo_numbers, random_rotate=False)
    olp = mol.obasis.compute_overlap()
    kin = mol.obasis.compute_kinetic()
    na = mol.obasis.compute_nuclear_attraction(mol.coordinates, mol.pseudo_numbers)
    terms = [
        RGridGroup(mol.obasis, grid, [
            RLibXCLDA('c_vwn'),
        ]),
    ]
    ham = REffHam(terms)
    return mol, olp, kin, na, ham


def test_cubic_interpolation_c_vwn_cs():
    mol, olp, kin, na, ham = setup_c_vwn_cs()
    check_interpolation(ham, olp, kin, na, [mol.orb_alpha])


def test_dot_hessian_c_vwn_cs():
    mol, _olp, _kin, _na, ham = setup_c_vwn_cs()
    check_dot_hessian(ham, mol.dm_alpha)


def test_dot_hessian_c_vwn_cs_polynomial():
    mol, olp, kin, na, ham = setup_c_vwn_cs()
    check_dot_hessian_polynomial(olp, kin+na, ham, [mol.orb_alpha], is_hf=False, extent=0.00005)


def test_dot_hessian_c_vwn_cs_cache():
    mol, _olp, _kin, _na, ham = setup_c_vwn_cs()
    check_dot_hessian_cache(ham, mol.dm_alpha)


def setup_o3lyp_cs():
    """Prepare datastructures for R-O3LYP (xc-functional-only) calculation on water."""
    fn_fchk = context.get_fn('test/water_hfs_321g.fchk')
    mol = IOData.from_file(fn_fchk)
    mol.dm_alpha = mol.orb_alpha.to_dm()

    grid = BeckeMolGrid(mol.coordinates, mol.numbers, mol.pseudo_numbers, random_rotate=False)
    olp = mol.obasis.compute_overlap()
    kin = mol.obasis.compute_kinetic()
    na = mol.obasis.compute_nuclear_attraction(mol.coordinates, mol.pseudo_numbers)
    libxc_term = RLibXCHybridGGA('xc_o3lyp')
    terms = [
        RGridGroup(mol.obasis, grid, [libxc_term]),
    ]
    ham = REffHam(terms)
    return mol, olp, kin, na, ham


def test_cubic_interpolation_o3lyp_cs():
    mol, olp, kin, na, ham = setup_o3lyp_cs()
    check_interpolation(ham, olp, kin, na, [mol.orb_alpha])


def test_dot_hessian_o3lyp_cs():
    mol, _olp, _kin, _na, ham = setup_o3lyp_cs()
    check_dot_hessian(ham, mol.dm_alpha)


def test_dot_hessian_o3lyp_cs_polynomial():
    raise SkipTest("We should use more robust tests for derivatives.")
    # mol, olp, kin, na, ham = setup_o3lyp_cs()
    # check_dot_hessian_polynomial(olp, kin+na, ham, [mol.orb_alpha], is_hf=False, extent=0.00001)


def test_dot_hessian_o3lyp_cs_cache():
    mol, _olp, _kin, _na, ham = setup_o3lyp_cs()
    check_dot_hessian_cache(ham, mol.dm_alpha)


def test_cubic_interpolation_x_tpss_cs():
    fn_fchk = context.get_fn('test/water_hfs_321g.fchk')
    mol = IOData.from_file(fn_fchk)

    grid = BeckeMolGrid(mol.coordinates, mol.numbers, mol.pseudo_numbers, random_rotate=False)
    olp = mol.obasis.compute_overlap()
    kin = mol.obasis.compute_kinetic()
    na = mol.obasis.compute_nuclear_attraction(mol.coordinates, mol.pseudo_numbers)
    terms = [
        RGridGroup(mol.obasis, grid, [RLibXCMGGA('x_tpss')]),
    ]
    ham = REffHam(terms)
    check_interpolation(ham, olp, kin, na, [mol.orb_alpha])


def test_cubic_interpolation_c_pbe_os():
    fn_fchk = context.get_fn('test/h3_pbe_321g.fchk')
    mol = IOData.from_file(fn_fchk)

    grid = BeckeMolGrid(mol.coordinates, mol.numbers, mol.pseudo_numbers, random_rotate=False)
    olp = mol.obasis.compute_overlap()
    kin = mol.obasis.compute_kinetic()
    na = mol.obasis.compute_nuclear_attraction(mol.coordinates, mol.pseudo_numbers)
    terms = [
        UGridGroup(mol.obasis, grid, [
            ULibXCGGA('c_pbe'),
        ]),
    ]
    ham = UEffHam(terms)
    check_interpolation(ham, olp, kin, na, [mol.orb_alpha, mol.orb_beta])


def test_cubic_interpolation_x_pbe_os():
    fn_fchk = context.get_fn('test/h3_pbe_321g.fchk')
    mol = IOData.from_file(fn_fchk)

    grid = BeckeMolGrid(mol.coordinates, mol.numbers, mol.pseudo_numbers, random_rotate=False)
    olp = mol.obasis.compute_overlap()
    kin = mol.obasis.compute_kinetic()
    na = mol.obasis.compute_nuclear_attraction(mol.coordinates, mol.pseudo_numbers)
    terms = [
        UGridGroup(mol.obasis, grid, [
            ULibXCGGA('x_pbe'),
        ]),
    ]
    ham = UEffHam(terms)
    check_interpolation(ham, olp, kin, na, [mol.orb_alpha, mol.orb_beta])


def setup_hfs_os():
    """Prepare datastructures for U_HFS (x-functional-only) calculation in H3 radical."""
    fn_fchk = context.get_fn('test/h3_hfs_321g.fchk')
    mol = IOData.from_file(fn_fchk)
    mol.dm_alpha = mol.orb_alpha.to_dm()
    mol.dm_beta = mol.orb_beta.to_dm()
    grid = BeckeMolGrid(mol.coordinates, mol.numbers, mol.pseudo_numbers, random_rotate=False)
    olp = mol.obasis.compute_overlap()
    kin = mol.obasis.compute_kinetic()
    na = mol.obasis.compute_nuclear_attraction(mol.coordinates, mol.pseudo_numbers)
    terms = [
        UGridGroup(mol.obasis, grid, [
            ULibXCLDA('x'),
        ]),
    ]
    ham = UEffHam(terms)
    return mol, olp, kin, na, ham


def test_cubic_interpolation_hfs_os():
    mol, olp, kin, na, ham = setup_hfs_os()
    check_interpolation(ham, olp, kin, na, [mol.orb_alpha, mol.orb_beta])


def test_cubic_interpolation_x_pbe_c_vwn_os():
    # mixing of LDA and GGA
    fn_fchk = context.get_fn('test/h3_hfs_321g.fchk')
    mol = IOData.from_file(fn_fchk)

    grid = BeckeMolGrid(mol.coordinates, mol.numbers, mol.pseudo_numbers, random_rotate=False)
    olp = mol.obasis.compute_overlap()
    kin = mol.obasis.compute_kinetic()
    na = mol.obasis.compute_nuclear_attraction(mol.coordinates, mol.pseudo_numbers)
    terms = [
        UGridGroup(mol.obasis, grid, [
            ULibXCGGA('x_pbe'),
            ULibXCLDA('c_vwn'),
        ]),
    ]
    ham = UEffHam(terms)
    check_interpolation(ham, olp, kin, na, [mol.orb_alpha, mol.orb_beta])


def test_cubic_interpolation_o3lyp_os():
    fn_fchk = context.get_fn('test/h3_hfs_321g.fchk')
    mol = IOData.from_file(fn_fchk)

    grid = BeckeMolGrid(mol.coordinates, mol.numbers, mol.pseudo_numbers, random_rotate=False)
    olp = mol.obasis.compute_overlap()
    kin = mol.obasis.compute_kinetic()
    na = mol.obasis.compute_nuclear_attraction(mol.coordinates, mol.pseudo_numbers)
    er = mol.obasis.compute_electron_repulsion()
    libxc_term = ULibXCHybridGGA('xc_o3lyp')
    terms = [
        UGridGroup(mol.obasis, grid, [libxc_term]),
        UExchangeTerm(er, 'x_hf', libxc_term.get_exx_fraction()),
    ]
    ham = UEffHam(terms)
    check_interpolation(ham, olp, kin, na, [mol.orb_alpha, mol.orb_beta])


def test_cubic_interpolation_x_tpss_os():
    # mixing of LDA and GGA
    fn_fchk = context.get_fn('test/h3_hfs_321g.fchk')
    mol = IOData.from_file(fn_fchk)

    grid = BeckeMolGrid(mol.coordinates, mol.numbers, mol.pseudo_numbers, random_rotate=False)
    olp = mol.obasis.compute_overlap()
    kin = mol.obasis.compute_kinetic()
    na = mol.obasis.compute_nuclear_attraction(mol.coordinates, mol.pseudo_numbers)
    terms = [
        UGridGroup(mol.obasis, grid, [ULibXCMGGA('x_tpss')]),
    ]
    ham = UEffHam(terms)
    check_interpolation(ham, olp, kin, na, [mol.orb_alpha, mol.orb_beta])


def test_functionals_present():
    t1 = RLibXCLDA('c_vwn')     # The VWN 5 functional
    assert t1._libxc_wrapper.key == 'lda_c_vwn'
    t2 = RLibXCLDA('c_vwn_4')   # The VWN 4 functional
    assert t2._libxc_wrapper.key == 'lda_c_vwn_4'
    t3 = RLibXCHybridGGA('xc_wb97x')
    assert t3._libxc_wrapper.key == 'hyb_gga_xc_wb97x'


ref_lda_x_1 = """\
@article{Dirac1930_376,
  title = {Note on Exchange Phenomena in the Thomas Atom},
  author = {P. A. M. Dirac},\n  journal = {Math. Proc. Cambridge Philos. Soc.},
  volume = {26},
  issue = {03},
  month = {7},
  pages = {376},
  year = {1930},
  issn = {1469-8064},
  doi = {10.1017/S0305004100016108},
  url = {http://journals.cambridge.org/article_S0305004100016108}
}"""


ref_lda_x_2 = """\
@article{Bloch1929_545,
  title = {Bemerkung zur Elektronentheorie des Ferromagnetismus und der elektrischen \
Leitf\\"{a}higkeit},
  author = {F. Bloch},
  journal = {Z. Phys.},
  volume = {57},
  number = {7-8},
  pages = {545},
  year = {1929},
  issn = {0044-3328},
  publisher = {Springer-Verlag},
  language = {German},
  doi = {10.1007/BF01340281},
  url = {http://link.springer.com/article/10.1007\\%2FBF01340281}
}"""


def test_info():
    t = RLibXCWrapper('lda_x')
    assert t.key == 'lda_x'
    assert t.name == "Slater exchange"
    assert isinstance(t.number, int)
    assert isinstance(t.kind, int)
    assert isinstance(t.family, int)
    assert t.refs == [
        ['P. A. M. Dirac, Math. Proc. Cambridge Philos. Soc. 26, 376 (1930)',
         '10.1017/S0305004100016108',
         ref_lda_x_1],
        ['F. Bloch, Z. Phys. 57, 545 (1929)',
         '10.1007/BF01340281',
         ref_lda_x_2]]


def test_info_nonexisting():
    with assert_raises(ValueError):
        RLibXCWrapper('lda_foobar')


def test_hyb_cam_exx_parameters():
    # xc_pbeh = The PBE0 functional
    t1 = RLibXCHybridGGA('xc_pbeh')
    assert t1.get_exx_fraction() == 0.25
    assert t1.get_cam_coeffs() == (0.0, 0.25, 0.0)
    t2 = ULibXCHybridGGA('xc_pbeh')
    assert t2.get_exx_fraction() == 0.25
    assert t2.get_cam_coeffs() == (0.0, 0.25, 0.0)
    t3 = RLibXCHybridGGA('xc_wb97x')
    assert t3.get_cam_coeffs() == (0.3, 1.0, -0.842294)
    assert t3.get_exx_fraction() == 1.0
    # xc_tpssh = The TPSS functional with exact exchange
    t4 = RLibXCHybridMGGA('xc_tpssh')
    assert t4.get_exx_fraction() == 0.1
    assert t4.get_cam_coeffs() == (0.0, 0.1, 0.0)
    t5 = ULibXCHybridMGGA('xc_tpssh')
    assert t5.get_exx_fraction() == 0.1
    assert t5.get_cam_coeffs() == (0.0, 0.1, 0.0)
