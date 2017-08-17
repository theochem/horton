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


import numpy as np

from .. proatomdb import ProAtomDB
from .common import get_fn, load_atoms_npz


def test_db_basics():
    records = load_atoms_npz(numbers=[8, 1], max_cation=1, max_anion=-1)
    padb = ProAtomDB(records)
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
    # assert padb.get_rgrid(8) is r1.rgrid
    assert padb.get_record(8, +1).ipot_energy is None
    assert padb.get_record(8, -1).ipot_energy == padb.get_record(8, 0).energy - padb.get_record(8, -1).energy
    assert padb.get_record(1, 0).ipot_energy == -padb.get_record(1, 0).energy


def test_db_basics_pseudo():
    records = load_atoms_npz(numbers=[8, 14], max_cation=2, max_anion=-2)
    padb = ProAtomDB(records)
    assert padb.get_numbers() == [8, 14]
    assert padb.get_charges(8) == [2, 1, 0, -1]
    assert padb.get_charges(8, safe=True) == [2, 1, 0, -1]
    assert padb.get_charges(14) == [2, 1, 0, -1]
    # assert padb.get_rgrid(8) is padb.get_record(8, -1).rgrid
    # assert padb.get_rgrid(8) is padb.get_record(8, 0).rgrid
    # assert padb.get_rgrid(8) is padb.get_record(8, 1).rgrid
    # assert padb.get_rgrid(8) is padb.get_record(8, 2).rgrid
    r1 = padb.get_record(8, -1)
    assert r1.safe
    assert abs(r1.energy - -72.5869876339) < 1e-8
    assert abs(r1.ipot_energy - (padb.get_record(8, 0).energy - r1.energy)) < 1e-5
    assert padb.get_record(8, +2).ipot_energy is None


def test_record_basics_pseudo():
    record = load_atoms_npz(numbers=[6], max_cation=0, max_anion=0)[0]
    assert record.number == 6
    assert record.charge == 0
    assert abs(record.energy - -37.7739982507) < 1e-10
    assert record.ipot_energy is None
    assert record.population == 6
    assert record.pseudo_number == 6
    assert record.pseudo_population == 6
    assert record.safe


def compare_padbs(padb1, padb2):
    assert padb1.size == padb2.size
    for number in padb1.get_numbers():
        for charge in padb1.get_charges(number):
            r1 = padb1.get_record(number, charge)
            r2 = padb2.get_record(number, charge)
            assert r1 == r2


def test_io_group():
    records = load_atoms_npz(numbers=[1, 6], max_cation=1, max_anion=-1)
    padb1 = ProAtomDB(records)
    assert padb1.size == 5
    keys = sorted(padb1._map.keys())
    assert keys == [(1, -1), (1, 0), (6, -1), (6, 0), (6, +1)]


def test_io_filename():
    records = load_atoms_npz(numbers=[1, 6], max_cation=1, max_anion=0)
    padb1 = ProAtomDB(records)
    keys = sorted(padb1._map.keys())
    assert keys == [(1, 0), (6, 0), (6, 1)]


def check_spline_record(spline, record):
    assert abs(spline.y - record.rho).max() < 1e-10
    assert abs(spline.dx - record.deriv).max() < 1e-10


def check_spline_pop(spline, pop):
    rtf = spline.rtransform
    check_pop = 4 * np.pi * np.sum(rtf.get_deriv() * rtf.get_radii()**2 * spline.y)
    assert abs(pop - check_pop) < 1e-2


def check_spline_mono_decr(spline):
    t = np.arange(0, spline.rtransform.npoint, 0.1)
    x = spline.rtransform.radius(t)
    y = spline(x)
    i = (abs(y) < 1e-10).nonzero()[0][0]
    y = y[:i]
    assert ((y[1:] - y[:-1]) / y[:-1]).min() < 1e-9


def test_get_spline():
    records = load_atoms_npz(numbers=[1, 6], max_cation=1, max_anion=-1)
    padb = ProAtomDB(records)

    spline = padb.get_spline(6)
    check_spline_pop(spline, 6.0)
    check_spline_record(spline, padb.get_record(6, 0))
    check_spline_mono_decr(spline)

    spline = padb.get_spline(6, -1)
    check_spline_pop(spline, 7.0)
    check_spline_record(spline, padb.get_record(6, -1))
    check_spline_mono_decr(spline)

    spline = padb.get_spline(6, {0: 0.5, -1: 0.5})
    check_spline_pop(spline, 6.5)
    check_spline_mono_decr(spline)

    spline = padb.get_spline(1, {0: 0.5})
    check_spline_pop(spline, 0.5)
    check_spline_mono_decr(spline)


def test_get_spline_pseudo():
    records = load_atoms_npz(numbers=[8, 14], max_cation=1, max_anion=-1)
    padb = ProAtomDB(records)

    spline = padb.get_spline(8)
    check_spline_pop(spline, 8.0)
    check_spline_record(spline, padb.get_record(8, 0))

    spline = padb.get_spline(8, -1)
    check_spline_pop(spline, 9.0)
    check_spline_record(spline, padb.get_record(8, -1))

    spline = padb.get_spline(8, {0: 0.5, -1: 0.5})
    check_spline_pop(spline, 8.5)

    spline = padb.get_spline(14)
    check_spline_pop(spline, 14.0)
    check_spline_record(spline, padb.get_record(14, 0))


def test_compact():
    records = load_atoms_npz(numbers=[8, 14], max_cation=1, max_anion=0)
    padb = ProAtomDB(records)
    padb.compact(0.1)
    assert padb.get_rgrid(8).size < 100
    assert padb.get_rgrid(14).size < 100


def test_normalize():
    records = load_atoms_npz(numbers=[8, 14], max_cation=1, max_anion=0)
    padb = ProAtomDB(records)
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
    records = load_atoms_npz(numbers=[8], max_cation=1, max_anion=0)
    padb = ProAtomDB(records)
    assert (padb.get_rho(8, {}) == 0.0).all()
