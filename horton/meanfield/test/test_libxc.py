# -*- coding: utf-8 -*-
# HORTON: Helpful Open-source Research TOol for N-fermion systems.
# Copyright (C) 2011-2017 The HORTON Development Team
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



from nose.plugins.skip import SkipTest
from nose.tools import assert_raises

from horton.grid import BeckeMolGrid
from .common import check_interpolation, \
    check_dot_hessian, check_dot_hessian_polynomial, check_dot_hessian_cache, load_olp, load_kin, \
    load_na, get_obasis, load_mdata, load_orbs_alpha, load_orbsa_dms, load_orbs_beta, \
    load_orbsb_dms, load_er
from .. import RGridGroup, RLibXCGGA, REffHam, RLibXCLDA, RLibXCHybridGGA, UGridGroup, ULibXCGGA, \
    RLibXCMGGA, UEffHam, ULibXCLDA, ULibXCHybridGGA, UExchangeTerm, ULibXCMGGA, RLibXCWrapper, \
    RLibXCHybridMGGA, ULibXCHybridMGGA


def setup_gga_cs(name):
    """Prepare data structures for R-GGA calculation in CO."""
    fname = 'co_pbe_sto3g_fchk'
    mdata = load_mdata(fname)
    dm_alpha = load_orbsa_dms(fname)
    orb_alpha = load_orbs_alpha(fname)

    grid = BeckeMolGrid(mdata['coordinates'], mdata['numbers'], mdata['pseudo_numbers'],
                        random_rotate=False)
    olp = load_olp(fname)
    kin = load_kin(fname)
    na = load_na(fname)
    terms = [
        RGridGroup(get_obasis(fname), grid, [
            RLibXCGGA(name),
        ]),
    ]
    ham = REffHam(terms)
    return dm_alpha, olp, kin, na, ham, orb_alpha


def test_cubic_interpolation_c_pbe_cs():
    dm_alpha, olp, kin, na, ham, orb_alpha = setup_gga_cs('c_pbe')
    check_interpolation(ham, olp, kin, na, [orb_alpha])


def test_dot_hessian_c_pbe_cs():
    dm_alpha, olp, kin, na, ham, orb_alpha = setup_gga_cs('c_pbe')
    check_dot_hessian(ham, dm_alpha)


def test_dot_hessian_c_pbe_cs_polynomial():
    dm_alpha, olp, kin, na, ham, orb_alpha = setup_gga_cs('c_pbe')
    check_dot_hessian_polynomial(olp, kin + na, ham, [orb_alpha], is_hf=False, extent=0.00005)


def test_dot_hessian_c_pbe_cs_cache():
    dm_alpha, olp, kin, na, ham, orb_alpha = setup_gga_cs('c_pbe')
    check_dot_hessian_cache(ham, dm_alpha)


def test_cubic_interpolation_x_pbe_cs():
    dm_alpha, olp, kin, na, ham, orb_alpha = setup_gga_cs('x_pbe')
    check_interpolation(ham, olp, kin, na, [orb_alpha])


def test_dot_hessian_x_pbe_cs():
    dm_alpha, olp, kin, na, ham, orb_alpha = setup_gga_cs('x_pbe')
    check_dot_hessian(ham, dm_alpha)


def test_dot_hessian_x_pbe_cs_polynomial():
    dm_alpha, olp, kin, na, ham, orb_alpha = setup_gga_cs('x_pbe')
    check_dot_hessian_polynomial(olp, kin + na, ham, [orb_alpha], is_hf=False, extent=0.00005)


def test_dot_hessian_x_pbe_cs_cache():
    dm_alpha, olp, kin, na, ham, orb_alpha = setup_gga_cs('x_pbe')
    check_dot_hessian_cache(ham, dm_alpha)


def setup_hfs_cs():
    """Prepare data structures for R-HFS (x-functional-only) calculation on CO."""
    fname = 'co_pbe_sto3g_fchk'
    mdata = load_mdata(fname)
    dm_alpha = load_orbsa_dms(fname)
    orb_alpha = load_orbs_alpha(fname)

    grid = BeckeMolGrid(mdata['coordinates'], mdata['numbers'], mdata['pseudo_numbers'],
                        'coarse', random_rotate=False)
    olp = load_olp(fname)
    kin = load_kin(fname)
    na = load_na(fname)
    terms = [
        RGridGroup(get_obasis(fname), grid, [
            RLibXCLDA('x'),
        ]),
    ]
    ham = REffHam(terms)

    return dm_alpha, olp, kin, na, ham, orb_alpha


def test_cubic_interpolation_hfs_cs():
    dm_alpha, olp, kin, na, ham, orb_alpha = setup_hfs_cs()
    check_interpolation(ham, olp, kin, na, [orb_alpha])


def test_dot_hessian_hfs_cs():
    dm_alpha, olp, kin, na, ham, orb_alpha = setup_hfs_cs()
    check_dot_hessian(ham, dm_alpha)


def test_dot_hessian_hfs_cs_polynomial():
    dm_alpha, olp, kin, na, ham, orb_alpha = setup_hfs_cs()
    check_dot_hessian_polynomial(olp, kin + na, ham, [orb_alpha], is_hf=False, extent=0.00005)


def test_dot_hessian_hfs_cs_cache():
    dm_alpha, olp, kin, na, ham, orb_alpha = setup_hfs_cs()
    check_dot_hessian_cache(ham, dm_alpha)


def setup_x_pbe_c_vwn_cs():
    """Setup data structure for mixed GGA+LDA calculation."""
    # mixing of GGA and LDA
    fname = 'water_hfs_321g_fchk'
    mdata = load_mdata(fname)
    dm_alpha = load_orbsa_dms(fname)
    orb_alpha = load_orbs_alpha(fname)

    grid = BeckeMolGrid(mdata['coordinates'], mdata['numbers'], mdata['pseudo_numbers'],
                        random_rotate=False)
    olp = load_olp(fname)
    kin = load_kin(fname)
    na = load_na(fname)
    terms = [
        RGridGroup(get_obasis(fname), grid, [
            RLibXCGGA('x_pbe'),
            RLibXCLDA('c_vwn'),
        ]),
    ]
    ham = REffHam(terms)
    return dm_alpha, olp, kin, na, ham, orb_alpha


def test_cubic_interpolation_x_pbe_c_vwn_cs():
    dm_alpha, olp, kin, na, ham, orb_alpha = setup_x_pbe_c_vwn_cs()
    check_interpolation(ham, olp, kin, na, [orb_alpha])


def test_dot_hessian_x_pbe_c_vwn_cs():
    dm_alpha, olp, kin, na, ham, orb_alpha = setup_x_pbe_c_vwn_cs()
    check_dot_hessian(ham, dm_alpha)


def test_dot_hessian_x_pbe_c_vwn_cs_polynomial():
    dm_alpha, olp, kin, na, ham, orb_alpha = setup_x_pbe_c_vwn_cs()
    check_dot_hessian_polynomial(olp, kin + na, ham, [orb_alpha], is_hf=False, extent=0.00005)


def test_dot_hessian_x_pbe_c_vwn_cs_cache():
    dm_alpha, olp, kin, na, ham, orb_alpha = setup_x_pbe_c_vwn_cs()
    check_dot_hessian_cache(ham, dm_alpha)


def setup_c_vwn_cs():
    """Prepare data structures for R-VWN (c-functional-only) calculation on water."""
    fname = 'water_hfs_321g_fchk'
    mdata = load_mdata(fname)
    dm_alpha = load_orbsa_dms(fname)
    orb_alpha = load_orbs_alpha(fname)

    grid = BeckeMolGrid(mdata['coordinates'], mdata['numbers'], mdata['pseudo_numbers'],
                        random_rotate=False)
    olp = load_olp(fname)
    kin = load_kin(fname)
    na = load_na(fname)
    terms = [
        RGridGroup(get_obasis(fname), grid, [
            RLibXCLDA('c_vwn'),
        ]),
    ]
    ham = REffHam(terms)
    return dm_alpha, olp, kin, na, ham, orb_alpha


def test_cubic_interpolation_c_vwn_cs():
    dm_alpha, olp, kin, na, ham, orb_alpha = setup_c_vwn_cs()
    check_interpolation(ham, olp, kin, na, [orb_alpha])


def test_dot_hessian_c_vwn_cs():
    dm_alpha, olp, kin, na, ham, orb_alpha = setup_c_vwn_cs()
    check_dot_hessian(ham, dm_alpha)


def test_dot_hessian_c_vwn_cs_polynomial():
    dm_alpha, olp, kin, na, ham, orb_alpha = setup_c_vwn_cs()
    check_dot_hessian_polynomial(olp, kin + na, ham, [orb_alpha], is_hf=False, extent=0.00005)


def test_dot_hessian_c_vwn_cs_cache():
    dm_alpha, olp, kin, na, ham, orb_alpha = setup_c_vwn_cs()
    check_dot_hessian_cache(ham, dm_alpha)


def setup_o3lyp_cs():
    """Prepare data structures for R-O3LYP (xc-functional-only) calculation on water."""
    fname = 'water_hfs_321g_fchk'
    mdata = load_mdata(fname)
    dm_alpha = load_orbsa_dms(fname)
    orb_alpha = load_orbs_alpha(fname)

    grid = BeckeMolGrid(mdata['coordinates'], mdata['numbers'], mdata['pseudo_numbers'],
                        random_rotate=False)
    olp = load_olp(fname)
    kin = load_kin(fname)
    na = load_na(fname)
    libxc_term = RLibXCHybridGGA('xc_o3lyp')
    terms = [
        RGridGroup(get_obasis(fname), grid, [libxc_term]),
    ]
    ham = REffHam(terms)
    return dm_alpha, olp, kin, na, ham, orb_alpha


def test_cubic_interpolation_o3lyp_cs():
    dm_alpha, olp, kin, na, ham, orb_alpha = setup_o3lyp_cs()
    check_interpolation(ham, olp, kin, na, [orb_alpha])


def test_dot_hessian_o3lyp_cs():
    dm_alpha, olp, kin, na, ham, orb_alpha = setup_o3lyp_cs()
    check_dot_hessian(ham, dm_alpha)


def test_dot_hessian_o3lyp_cs_polynomial():
    raise SkipTest("We should use more robust tests for derivatives.")
    # dm_alpha, olp, kin, na, ham, orb_alpha = setup_o3lyp_cs()
    # check_dot_hessian_polynomial(olp, kin+na, ham, [orb_alpha], is_hf=False, extent=0.00001)


def test_dot_hessian_o3lyp_cs_cache():
    dm_alpha, olp, kin, na, ham, orb_alpha = setup_o3lyp_cs()
    check_dot_hessian_cache(ham, dm_alpha)


def test_cubic_interpolation_x_tpss_cs():
    fname = 'water_hfs_321g_fchk'
    mdata = load_mdata(fname)

    grid = BeckeMolGrid(mdata['coordinates'], mdata['numbers'], mdata['pseudo_numbers'],
                        random_rotate=False)
    olp = load_olp(fname)
    kin = load_kin(fname)
    na = load_na(fname)
    terms = [
        RGridGroup(get_obasis(fname), grid, [RLibXCMGGA('x_tpss')]),
    ]
    ham = REffHam(terms)
    check_interpolation(ham, olp, kin, na, [load_orbs_alpha(fname)])


def test_cubic_interpolation_c_pbe_os():
    fname = 'h3_pbe_321g_fchk'
    mdata = load_mdata(fname)

    grid = BeckeMolGrid(mdata['coordinates'], mdata['numbers'], mdata['pseudo_numbers'],
                        random_rotate=False)
    olp = load_olp(fname)
    kin = load_kin(fname)
    na = load_na(fname)
    terms = [
        UGridGroup(get_obasis(fname), grid, [
            ULibXCGGA('c_pbe'),
        ]),
    ]
    ham = UEffHam(terms)
    check_interpolation(ham, olp, kin, na, [load_orbs_alpha(fname), load_orbs_beta(fname)])


def test_cubic_interpolation_x_pbe_os():
    fname = 'h3_pbe_321g_fchk'
    mdata = load_mdata(fname)

    grid = BeckeMolGrid(mdata['coordinates'], mdata['numbers'], mdata['pseudo_numbers'],
                        random_rotate=False)
    olp = load_olp(fname)
    kin = load_kin(fname)
    na = load_na(fname)
    terms = [
        UGridGroup(get_obasis(fname), grid, [
            ULibXCGGA('x_pbe'),
        ]),
    ]
    ham = UEffHam(terms)
    check_interpolation(ham, olp, kin, na, [load_orbs_alpha(fname), load_orbs_beta(fname)])


def setup_hfs_os():
    """Prepare data structures for U_HFS (x-functional-only) calculation in H3 radical."""
    fname = 'h3_hfs_321g_fchk'
    mdata = load_mdata(fname)
    dm_alpha = load_orbsa_dms(fname)
    dm_beta = load_orbsb_dms(fname)
    orb_alpha = load_orbs_alpha(fname)
    orb_beta = load_orbs_beta(fname)
    grid = BeckeMolGrid(mdata['coordinates'], mdata['numbers'], mdata['pseudo_numbers'],
                        random_rotate=False)
    olp = load_olp(fname)
    kin = load_kin(fname)
    na = load_na(fname)
    terms = [
        UGridGroup(get_obasis(fname), grid, [
            ULibXCLDA('x'),
        ]),
    ]
    ham = UEffHam(terms)
    return dm_alpha, dm_beta, olp, kin, na, ham, orb_alpha, orb_beta


def test_cubic_interpolation_hfs_os():
    dm_alpha, dm_beta, olp, kin, na, ham, orb_alpha, orb_beta = setup_hfs_os()
    check_interpolation(ham, olp, kin, na, [orb_alpha, orb_beta])


def test_cubic_interpolation_x_pbe_c_vwn_os():
    # mixing of LDA and GGA
    fname = 'h3_hfs_321g_fchk'
    mdata = load_mdata(fname)

    grid = BeckeMolGrid(mdata['coordinates'], mdata['numbers'], mdata['pseudo_numbers'],
                        random_rotate=False)
    olp = load_olp(fname)
    kin = load_kin(fname)
    na = load_na(fname)
    terms = [
        UGridGroup(get_obasis(fname), grid, [
            ULibXCGGA('x_pbe'),
            ULibXCLDA('c_vwn'),
        ]),
    ]
    ham = UEffHam(terms)
    check_interpolation(ham, olp, kin, na, [load_orbs_alpha(fname), load_orbs_beta(fname)])


def test_cubic_interpolation_o3lyp_os():
    fname = 'h3_hfs_321g_fchk'
    mdata = load_mdata(fname)

    grid = BeckeMolGrid(mdata['coordinates'], mdata['numbers'], mdata['pseudo_numbers'],
                        random_rotate=False)
    olp = load_olp(fname)
    kin = load_kin(fname)
    na = load_na(fname)
    er = load_er(fname)
    libxc_term = ULibXCHybridGGA('xc_o3lyp')
    terms = [
        UGridGroup(get_obasis(fname), grid, [libxc_term]),
        UExchangeTerm(er, 'x_hf', libxc_term.get_exx_fraction()),
    ]
    ham = UEffHam(terms)
    check_interpolation(ham, olp, kin, na, [load_orbs_alpha(fname), load_orbs_beta(fname)])


def test_cubic_interpolation_x_tpss_os():
    # mixing of LDA and GGA
    fname = 'h3_hfs_321g_fchk'
    mdata = load_mdata(fname)

    grid = BeckeMolGrid(mdata['coordinates'], mdata['numbers'], mdata['pseudo_numbers'],
                        random_rotate=False)
    olp = load_olp(fname)
    kin = load_kin(fname)
    na = load_na(fname)
    terms = [
        UGridGroup(get_obasis(fname), grid, [ULibXCMGGA('x_tpss')]),
    ]
    ham = UEffHam(terms)
    check_interpolation(ham, olp, kin, na, [load_orbs_alpha(fname), load_orbs_beta(fname)])


def test_functionals_present():
    t1 = RLibXCLDA('c_vwn')  # The VWN 5 functional
    assert t1._libxc_wrapper.key == 'lda_c_vwn'
    t2 = RLibXCLDA('c_vwn_4')  # The VWN 4 functional
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
  numpages = {10},
  year = {1930},
  issn = {1469-8064},
  doi = {10.1017/S0305004100016108},
  URL = {http://journals.cambridge.org/article_S0305004100016108}
}"""

ref_lda_x_2 = """\
@article{Bloch1929_545,
  title = {Bemerkung zur Elektronentheorie des Ferromagnetismus und der elektrischen \
Leitf\xc3\xa4higkeit},
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
