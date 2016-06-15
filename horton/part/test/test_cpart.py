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

from horton import *  # pylint: disable=wildcard-import,unused-wildcard-import
from horton.part.test.common import check_names, check_proatom_splines, \
    get_fake_co, get_fake_pseudo_oo
from horton.scripts.cpart import cpart_schemes


def check_jbw_coarse(local):
    # This test is not supposed to generate meaningful numbers. The cube data
    # is too coarse and the reference atoms may have little similarities with
    # the DFT density.

    # Load the cube file
    fn_cube = context.get_fn('test/jbw_coarse_aedens.cube')
    mol = IOData.from_file(fn_cube)
    moldens = mol.cube_data
    ugrid = mol.grid

    # Load some pro-atoms
    proatomdb = ProAtomDB.from_refatoms(numbers=[8,14], max_kation=0, max_anion=0)

    # Run the partitioning
    cpart = HirshfeldCPart(mol.coordinates, mol.numbers, mol.pseudo_numbers,
                           ugrid, moldens, proatomdb, local=local,
                           wcor_numbers=range(119))

    # Do some testing
    assert cpart.local == local
    if local:
        names = cpart.do_all()
        check_names(names, cpart)
    else:
        cpart.do_charges()
        wcor = cpart.get_wcor()
        assert abs(cpart['pseudo_populations'].sum() - ugrid.integrate(wcor, moldens)) < 1e-10


@attr('slow')
def test_hirshfeld_jbw_coarse_local():
    check_jbw_coarse(True)


def test_hirshfeld_jbw_coarse_global():
    check_jbw_coarse(False)


def check_fake(scheme, pseudo, dowcor, local, absmean, **kwargs):
    if pseudo:
        coordinates, numbers, pseudo_numbers, ugrid, moldens, proatomdb = get_fake_pseudo_oo()
    else:
        coordinates, numbers, ugrid, moldens, proatomdb = get_fake_co()
        pseudo_numbers = numbers.astype(float)

    if dowcor:
        wcor_numbers = range(119)
    else:
        wcor_numbers = []

    CPartClass = cpart_schemes[scheme]
    cpart = CPartClass(coordinates, numbers, pseudo_numbers, ugrid,
                       moldens, proatomdb, local=local,
                       wcor_numbers=wcor_numbers, **kwargs)

    assert cpart.local == local
    cpart.do_charges()
    charges = cpart['charges']
    #print abs(charges.sum()), abs(charges).mean(), abs(abs(charges).mean() - absmean)
    assert abs(charges.sum()) < 1e-2
    assert abs(abs(charges).mean() - absmean) < 1e-3

    if kwargs.get('greedy', False):
        # In case of a greedy algorithm, one should compare the proatoms
        # in the partitioning object with proatoms directly evaluate with the
        # proatom splines
        check_proatom_splines(cpart)


@attr('slow')
def test_hirshfeld_fake_local():
    check_fake('h', pseudo=False, dowcor=True, local=True, absmean=0.112)


@attr('slow')
def test_hirshfeld_fake_global():
    check_fake('h', pseudo=False, dowcor=True, local=False, absmean=0.112)


def test_hirshfeld_fake_pseudo_local():
    check_fake('h', pseudo=True, dowcor=True, local=True, absmean=0.213)


def test_hirshfeld_fake_pseudo_global():
    check_fake('h', pseudo=True, dowcor=True, local=False, absmean=0.213)


def test_hirshfeld_fake_pseudo_nowcor_local():
    check_fake('h', pseudo=True, dowcor=True, local=True, absmean=0.213)


def test_hirshfeld_fake_pseudo_nowcor_global():
    check_fake('h', pseudo=True, dowcor=True, local=False, absmean=0.213)



@attr('slow')
def test_hirshfeld_i_fake_local():
    check_fake('hi', pseudo=False, dowcor=True, local=True, absmean=0.428, threshold=1e-5)


@attr('slow')
def test_hirshfeld_i_fake_global():
    check_fake('hi', pseudo=False, dowcor=True, local=False, absmean=0.428, threshold=1e-5)


@attr('slow')
def test_hirshfeld_i_fake_pseudo_local():
    check_fake('hi', pseudo=True, dowcor=True, local=True, absmean=0.400, threshold=1e-4)


def test_hirshfeld_i_fake_pseudo_global():
    check_fake('hi', pseudo=True, dowcor=True, local=False, absmean=0.400, threshold=1e-4)


@attr('slow')
def test_hirshfeld_i_fake_pseudo_nowcor_local():
    check_fake('hi', pseudo=True, dowcor=True, local=True, absmean=0.400, threshold=1e-4)


def test_hirshfeld_i_fake_pseudo_nowcor_global():
    check_fake('hi', pseudo=True, dowcor=True, local=False, absmean=0.400, threshold=1e-4)


@attr('slow')
def test_hirshfeld_i_fake_local_greedy():
    check_fake('hi', pseudo=False, dowcor=True, local=True, absmean=0.428, threshold=1e-5, greedy=True)


@attr('slow')
def test_hirshfeld_i_fake_global_greedy():
    check_fake('hi', pseudo=False, dowcor=True, local=False, absmean=0.428, threshold=1e-5, greedy=True)


@attr('slow')
def test_hirshfeld_i_fake_pseudo_local_greedy():
    check_fake('hi', pseudo=True, dowcor=True, local=True, absmean=0.400, threshold=1e-4, greedy=True)


def test_hirshfeld_i_fake_pseudo_global_greedy():
    check_fake('hi', pseudo=True, dowcor=True, local=False, absmean=0.400, threshold=1e-4, greedy=True)


@attr('slow')
def test_hirshfeld_i_fake_pseudo_nowcor_local_greedy():
    check_fake('hi', pseudo=True, dowcor=True, local=True, absmean=0.400, threshold=1e-4, greedy=True)


def test_hirshfeld_i_fake_pseudo_nowcor_global_greedy():
    check_fake('hi', pseudo=True, dowcor=True, local=False, absmean=0.400, threshold=1e-4, greedy=True)




@attr('slow')
def test_hirshfeld_e_fake_local():
    check_fake('he', pseudo=False, dowcor=True, local=True, absmean=0.323, threshold=1e-4)


@attr('slow')
def test_hirshfeld_e_fake_global():
    check_fake('he', pseudo=False, dowcor=True, local=False, absmean=0.373, threshold=1e-4)


@attr('slow')
def test_hirshfeld_e_fake_pseudo_local():
    check_fake('he', pseudo=True, dowcor=True, local=True, absmean=0.396, threshold=1e-4)


@attr('slow')
def test_hirshfeld_e_fake_pseudo_global():
    check_fake('he', pseudo=True, dowcor=True, local=False, absmean=0.396, threshold=1e-4)


@attr('slow')
def test_hirshfeld_e_fake_pseudo_nowcor_local():
    check_fake('he', pseudo=True, dowcor=True, local=True, absmean=0.396, threshold=1e-4)


@attr('slow')
def test_hirshfeld_e_fake_pseudo_nowcor_global():
    check_fake('he', pseudo=True, dowcor=True, local=False, absmean=0.396, threshold=1e-4)


@attr('slow')
def test_hirshfeld_e_fake_local_greedy():
    check_fake('he', pseudo=False, dowcor=True, local=True, absmean=0.323, threshold=1e-4, greedy=True)


@attr('slow')
def test_hirshfeld_e_fake_global_greedy():
    check_fake('he', pseudo=False, dowcor=True, local=False, absmean=0.374, threshold=1e-4, greedy=True)


@attr('slow')
def test_hirshfeld_e_fake_pseudo_local_greedy():
    check_fake('he', pseudo=True, dowcor=True, local=True, absmean=0.396, threshold=1e-4, greedy=True)


@attr('slow')
def test_hirshfeld_e_fake_pseudo_global_greedy():
    check_fake('he', pseudo=True, dowcor=True, local=False, absmean=0.396, threshold=1e-4, greedy=True)


@attr('slow')
def test_hirshfeld_e_fake_pseudo_nowcor_local_greedy():
    check_fake('he', pseudo=True, dowcor=True, local=True, absmean=0.396, threshold=1e-4, greedy=True)


@attr('slow')
def test_hirshfeld_e_fake_pseudo_nowcor_global_greedy():
    check_fake('he', pseudo=True, dowcor=True, local=False, absmean=0.396, threshold=1e-4, greedy=True)
