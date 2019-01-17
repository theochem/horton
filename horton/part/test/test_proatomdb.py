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


import numpy as np, h5py as h5

from horton import *  # pylint: disable=wildcard-import,unused-wildcard-import
from horton.part.test.common import get_proatomdb_cp2k
from horton.test.common import tmpdir


def test_db_basics():
    padb = ProAtomDB.from_refatoms(numbers=[8, 1], max_cation=1, max_anion=1)
    assert padb.get_numbers() == [1, 8]
    assert padb.get_charges(8) == [1, 0, -1]
    assert padb.get_charges(1) == [0, -1]
    r1 = padb.get_record(8, -1)
    assert r1.number == 8
    assert r1.charge == -1
    assert abs(r1.energy - -72.587) < 1e-3
    assert r1.ipot_energy == padb.get_record(8, 0).energy - r1.energy
    assert r1.population == 9
    assert r1.pseudo_number == 8
    assert r1.pseudo_population == 9
    assert r1.safe
    assert r1.rgrid.size == 59
    r2 = padb.get_record(8, -1)
    r3 = padb.get_record(8, 0)
    assert r1 == r2
    assert r1 != r3
    assert padb.get_rgrid(8) is r1.rgrid
    assert padb.get_record(8, +1).ipot_energy is None
    assert padb.get_record(8, -1).ipot_energy == padb.get_record(8, 0).energy - padb.get_record(8, -1).energy
    assert padb.get_record(1, 0).ipot_energy == -padb.get_record(1, 0).energy


def test_db_basics_pseudo():
    padb = get_proatomdb_cp2k()
    assert padb.get_numbers() == [8, 14]
    assert padb.get_charges(8) == [2, 1, 0, -1, -2]
    assert padb.get_charges(8, safe=True) == [2, 1, 0, -1]
    assert padb.get_charges(14) == [0]
    assert not padb.get_record(8, -2).safe
    assert padb.get_rgrid(8) is padb.get_record(8, -2).rgrid
    assert padb.get_rgrid(8) is padb.get_record(8, -1).rgrid
    assert padb.get_rgrid(8) is padb.get_record(8, 0).rgrid
    assert padb.get_rgrid(8) is padb.get_record(8, 1).rgrid
    assert padb.get_rgrid(8) is padb.get_record(8, 2).rgrid
    r1 = padb.get_record(8, -1)
    assert r1.safe
    assert abs(r1.energy - -15.866511882272) < 1e-8
    assert abs(r1.ipot_energy - (padb.get_record(8, 0).energy - r1.energy)) < 1e-5
    r2 = padb.get_record(8, -2)
    assert not r2.safe
    assert abs(r2.energy - -15.464982778766) < 1e-8
    assert abs(r2.ipot_energy - (r1.energy - r2.energy)) < 1e-5
    assert padb.get_record(8, +2).ipot_energy is None


def test_record_basics_pseudo():
    fn_out = context.get_fn('test/atom_si.cp2k.out')
    mol = IOData.from_file(fn_out)
    r = ProAtomRecord.from_iodata(mol)
    assert r.number == 14
    assert r.charge == 0
    assert abs(r.energy - -3.761587698067) < 1e-10
    assert r.ipot_energy is None
    assert r.population == 14
    assert r.pseudo_number == 4
    assert r.pseudo_population == 4
    assert r.safe


def compare_padbs(padb1, padb2):
    assert padb1.size == padb2.size
    for number in padb1.get_numbers():
        for charge in padb1.get_charges(number):
            r1 = padb1.get_record(number, charge)
            r2 = padb2.get_record(number, charge)
            assert r1 == r2


def test_io_group():
    padb1 = ProAtomDB.from_refatoms(numbers=[1, 6], max_cation=1, max_anion=1)
    assert padb1.size == 5
    keys = sorted(padb1._map.keys())
    assert keys == [(1, -1), (1, 0), (6, -1), (6, 0), (6, +1)]

    with h5.File('horton.dpart.test.test_proatomdb.test_io_group', driver='core', backing_store=False) as f:
        padb1.to_file(f)
        padb2 = ProAtomDB.from_file(f)
        compare_padbs(padb1, padb2)


def test_io_filename():
    padb1 = ProAtomDB.from_refatoms(numbers=[1, 6], max_cation=1, max_anion=0)
    keys = sorted(padb1._map.keys())
    assert keys == [(1, 0), (6, 0), (6, 1)]

    with tmpdir('horton.dpart.test.test_proatomdb.test_io_filename') as dn:
        filename = '%s/test.h5' % dn
        padb1.to_file(filename)
        padb2 = ProAtomDB.from_file(filename)
        compare_padbs(padb1, padb2)


def test_compute_radii():
    rgrid = RadialGrid(ExpRTransform(1e-3, 1e1, 100))
    padb = ProAtomDB.from_refatoms([1, 6], 0, 0, (rgrid, 110))
    record = padb.get_record(6, 0)
    indexes, radii = record.compute_radii([2.0, 5.9, 5.999])
    assert (indexes == [69, 90, 100]).all()
    assert abs(radii - np.array([0.600577, 4.168655, 10.0])).max() < 1e-5


def test_moments():
    padb = get_proatomdb_cp2k()
    record0 = padb.get_record(8, 0)
    record1 = padb.get_record(8, 1)
    m0 = record0.get_moment(3)
    m1 = record1.get_moment(3)
    assert m0 > m1
    assert abs(m0-21.84) < 1e-2
    assert abs(m1-12.17) < 1e-2


def check_spline_record(spline, record):
    assert abs(spline.y - record.rho).max() < 1e-10
    assert abs(spline.dx - record.deriv).max() < 1e-10


def check_spline_pop(spline, pop):
    rtf = spline.rtransform
    int1d = spline.rtransform.get_default_int1d()
    check_pop = 4*np.pi*dot_multi(
        rtf.get_deriv(),
        rtf.get_radii()**2,
        spline.y,
        int1d.get_weights(rtf.npoint),
    )
    assert abs(pop - check_pop) < 1e-2


def check_spline_mono_decr(spline):
    t = np.arange(0, spline.rtransform.npoint, 0.1)
    x = spline.rtransform.radius(t)
    y = spline(x)
    i = (abs(y) < 1e-10).nonzero()[0][0]
    y = y[:i]
    assert ((y[1:] - y[:-1])/y[:-1]).min() < 1e-9


def test_get_spline():
    padb = ProAtomDB.from_refatoms(numbers=[1, 6], max_cation=1, max_anion=1)

    spline = padb.get_spline(6)
    check_spline_pop(spline, 6.0)
    check_spline_record(spline, padb.get_record(6, 0))
    check_spline_mono_decr(spline)

    spline = padb.get_spline(6, -1)
    check_spline_pop(spline, 7.0)
    check_spline_record(spline, padb.get_record(6, -1))
    check_spline_mono_decr(spline)

    spline = padb.get_spline(6, {0:0.5, -1:0.5})
    check_spline_pop(spline, 6.5)
    check_spline_mono_decr(spline)

    spline = padb.get_spline(1, {0:0.5})
    check_spline_pop(spline, 0.5)
    check_spline_mono_decr(spline)


def test_get_spline_pseudo():
    padb = get_proatomdb_cp2k()

    spline = padb.get_spline(8)
    check_spline_pop(spline, 6.0)
    check_spline_record(spline, padb.get_record(8, 0))

    spline = padb.get_spline(8, -1)
    check_spline_pop(spline, 7.0)
    check_spline_record(spline, padb.get_record(8, -1))

    spline = padb.get_spline(8, {0:0.5, -1:0.5})
    check_spline_pop(spline, 6.5)

    spline = padb.get_spline(14)
    check_spline_pop(spline, 4.0)
    check_spline_record(spline, padb.get_record(14, 0))


def test_compact():
    padb = get_proatomdb_cp2k()
    padb.compact(0.1)
    assert padb.get_rgrid(8).size < 100
    assert padb.get_rgrid(14).size < 100


def test_normalize():
    padb = get_proatomdb_cp2k()
    padb.compact(0.1)
    padb.normalize()
    for number in padb.get_numbers():
        rgrid = padb.get_rgrid(number)
        for charge in padb.get_charges(number):
            r = padb.get_record(number, charge)
            nel = rgrid.integrate(r.rho)
            nel_integer = r.pseudo_number - charge
            assert abs(nel - nel_integer) < 1e-10


def test_empty_proatom():
    padb = get_proatomdb_cp2k()
    assert (padb.get_rho(8, {}) == 0.0).all()


def test_io_atdens():
    padb = ProAtomDB.from_file(context.get_fn('test/pro.atdens'))
    assert padb.get_numbers() == [16]
    assert padb.get_charges(16) == [3, 2]
    r = padb.get_record(16, 3)
    assert abs(r.rho[0] - 0.2628105459E+04) < 1e-5
    assert abs(r.rho[-1] - 0.1998952826E-16) < 1e-5
    s = padb.get_spline(16, 3)
    assert abs(s(np.array([0.0])) - 2661.68659449) < 1e-5
    radii = r.rgrid.rtransform.get_radii()
    assert radii[0] == 0.5216488380E-03
    assert abs(radii[-1] - 20) < 1e-14
    assert abs(radii[1] - 0.5442350204E-03) < 1e-8
    assert abs(r.rgrid.integrate(r.rho) - 13) < 1e-3
    # check the basics of the get_rho method (charge)
    rho1 = padb.get_rho(16, 3)
    rho2, deriv = padb.get_rho(16, 3, do_deriv=True)
    assert (rho1 == rho2).all()
    assert deriv is None
    # check the basics of the get_rho method (dict)
    rho1 = padb.get_rho(16, {3:1})
    rho2, deriv = padb.get_rho(16, {3:1}, do_deriv=True)
    assert (rho1 == rho2).all()
    assert deriv is None
