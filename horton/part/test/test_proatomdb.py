# -*- coding: utf-8 -*-
# Horton is a Density Functional Theory program.
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


import numpy as np, h5py as h5
import tempfile, os

from horton import *
from horton.part.test.common import get_proatomdb_ref, get_proatomdb_cp2k


def test_db_basics():
    padb = get_proatomdb_ref([8, 1], 1, 1)
    assert padb.get_numbers() == [1, 8]
    assert padb.get_charges(8) == [1, 0, -1]
    assert padb.get_charges(1) == [0, -1]
    r1 = padb.get_record(8, -1)
    assert r1.number == 8
    assert r1.charge == -1
    assert abs(r1.energy - -72.587) < 1e-3
    assert r1.population == 9
    assert r1.pseudo_number == 8
    assert r1.pseudo_population == 9
    assert r1.safe
    assert r1.rgrid.size == 100
    r2 = padb.get_record(8, -1)
    r3 = padb.get_record(8, 0)
    assert r1 == r2
    assert r1 != r3
    assert padb.get_rgrid(8) is r1.rgrid


def test_db_basics_pseudo():
    padb = get_proatomdb_cp2k()
    assert padb.get_numbers() == [8, 14]
    assert padb.get_charges(8) == [2, 1, 0, -1, -2]
    assert padb.get_charges(8, safe=True) == [2, 1, 0, -1]
    assert padb.get_charges(14) == [0]
    assert padb.get_record(8, -1).safe
    assert not padb.get_record(8, -2).safe
    assert padb.get_rgrid(8) is padb.get_record(8, 1).rgrid


def test_record_basics_pseudo():
    fn_out = context.get_fn('test/atom_si.cp2k.out')
    sys = System.from_file(fn_out)

    rtf = ExpRTransform(1e-3, 1e1, 100)
    rgrid = RadialIntGrid(rtf)
    atgrid = AtomicGrid(0, np.zeros(3, float), (rgrid, 110), random_rotate=False)

    r = ProAtomRecord.from_system(sys, atgrid)
    assert r.number == 14
    assert r.charge == 0
    assert abs(r.energy - -3.761587698067) < 1e-10
    assert r.population == 14
    assert r.pseudo_number == 4
    assert r.pseudo_population == 4
    assert r.safe
    assert r.rgrid.rtransform is rtf


def compare_padbs(padb1, padb2):
    assert padb1.size == padb2.size
    for number in padb1.get_numbers():
        for charge in padb1.get_charges(number):
            r1 = padb1.get_record(number, charge)
            r2 = padb2.get_record(number, charge)
            assert r1 == r2


def test_io_group():
    padb1 = get_proatomdb_ref([1, 6], 1, 1)
    assert padb1.size == 5
    keys = sorted(padb1._map.keys())
    assert keys == [(1, -1), (1, 0), (6, -1), (6, 0), (6, +1)]

    with h5.File('horton.dpart.test.test_proatomdb.test_io_group', driver='core', backing_store=False) as f:
        padb1.to_file(f)
        padb2 = ProAtomDB.from_file(f)
        compare_padbs(padb1, padb2)


def test_io_filename():
    padb1 = get_proatomdb_ref([1, 6], 1, 0)
    keys = sorted(padb1._map.keys())
    assert keys == [(1, 0), (6, 0), (6, 1)]

    tmpdir = tempfile.mkdtemp('horton.dpart.test.test_proatomdb.test_io_filename')
    filename = '%s/test.h5' % tmpdir
    try:
        padb1.to_file(filename)
        padb2 = ProAtomDB.from_file(filename)
        compare_padbs(padb1, padb2)
    finally:
        if os.path.isfile(filename):
            os.remove(filename)
        os.rmdir(tmpdir)


def test_compute_radii():
    padb = get_proatomdb_ref([6], 0, 0)
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
    assert abs(spline.copy_y() - record.rho).max() < 1e-10


def check_spline_pop(spline, pop):
    int1d = SimpsonIntegrator1D()
    rtf = spline.rtransform
    check_pop = 4*np.pi*dot_multi(
        rtf.get_volume_elements(),
        rtf.get_radii()**2,
        spline.copy_y(),
        int1d.get_weights(rtf.npoint),
    )
    assert abs(pop - check_pop) < 1e-2


def test_get_spline():
    padb = get_proatomdb_ref([1, 6], 1, 1)

    spline = padb.get_spline(6)
    check_spline_pop(spline, 6.0)
    check_spline_record(spline, padb.get_record(6, 0))

    spline = padb.get_spline(6, -1)
    check_spline_pop(spline, 7.0)
    check_spline_record(spline, padb.get_record(6, -1))

    spline = padb.get_spline(6, {0:0.5, -1:0.5})
    check_spline_pop(spline, 6.5)

    spline = padb.get_spline(1, {0:0.5})
    check_spline_pop(spline, 0.5)


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
