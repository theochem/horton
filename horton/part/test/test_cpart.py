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


from horton import *
from horton.part.test.common import check_names, check_proatom_splines, \
    get_fake_co, get_fake_pseudo_oo


def check_jbw_coarse(local):
    # This test is not supposed to generate meaningful numbers. The cube data
    # is too coarse and the reference atoms may have little similarities with
    # the DFT density.

    # Load the cube file
    fn_cube = context.get_fn('test/jbw_coarse_aedens.cube')
    sys = System.from_file(fn_cube)
    mol_dens = sys.extra['cube_data']
    ugrid = sys.grid

    # Load some pro-atoms
    proatomdb = ProAtomDB.from_refatoms(numbers=[8,14], max_kation=0, max_anion=0)

    # Run the partitioning
    cpart = HirshfeldCPart(sys, ugrid, local, mol_dens, proatomdb, range(119))

    # Do some testing
    if local:
        names = cpart.do_all()
        check_names(names, cpart)
    else:
        cpart.do_charges()
        wcor = cpart.get_wcor()
        assert abs(cpart['populations'].sum() - ugrid.integrate(wcor, mol_dens)) < 1e-10


def test_hirshfeld_jbw_coarse_local():
    check_jbw_coarse(True)


def test_hirshfeld_jbw_coarse_global():
    check_jbw_coarse(False)


def check_fake(scheme, pseudo, dowcor, local, absmean, **kwargs):
    if pseudo:
        sys, ugrid, mol_dens, proatomdb = get_fake_pseudo_oo()
    else:
        sys, ugrid, mol_dens, proatomdb = get_fake_co()

    if dowcor:
        wcor_numbers = range(119)
    else:
        wcor_numbers = []

    CPartClass = cpart_schemes[scheme]
    cpart = CPartClass(sys, ugrid, local, mol_dens, proatomdb, wcor_numbers, **kwargs)

    cpart.do_charges()
    charges = cpart['charges']
    print abs(charges.sum()), abs(charges).mean(), abs(abs(charges).mean() - absmean)
    assert abs(charges.sum()) < 1e-2
    assert abs(abs(charges).mean() - absmean) < 1e-3

    if kwargs.get('greedy', False):
        # In case of a greedy algorithm, one should compare the proatoms
        # in the partitioning object with proatoms directly evaluate with the
        # proatom splines
        check_proatom_splines(cpart)


def test_hirshfeld_fake_local():
    check_fake('h', pseudo=False, dowcor=True, local=True, absmean=0.112)


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



def test_hirshfeld_i_fake_local():
    check_fake('hi', pseudo=False, dowcor=True, local=True, absmean=0.428, threshold=1e-5)


def test_hirshfeld_i_fake_global():
    check_fake('hi', pseudo=False, dowcor=True, local=False, absmean=0.428, threshold=1e-5)


def test_hirshfeld_i_fake_pseudo_local():
    check_fake('hi', pseudo=True, dowcor=True, local=True, absmean=0.400, threshold=1e-4)


def test_hirshfeld_i_fake_pseudo_global():
    check_fake('hi', pseudo=True, dowcor=True, local=False, absmean=0.400, threshold=1e-4)


def test_hirshfeld_i_fake_pseudo_nowcor_local():
    check_fake('hi', pseudo=True, dowcor=True, local=True, absmean=0.400, threshold=1e-4)


def test_hirshfeld_i_fake_pseudo_nowcor_global():
    check_fake('hi', pseudo=True, dowcor=True, local=False, absmean=0.400, threshold=1e-4)


def test_hirshfeld_i_fake_local_greedy():
    check_fake('hi', pseudo=False, dowcor=True, local=True, absmean=0.428, threshold=1e-5, greedy=True)


def test_hirshfeld_i_fake_global_greedy():
    check_fake('hi', pseudo=False, dowcor=True, local=False, absmean=0.428, threshold=1e-5, greedy=True)


def test_hirshfeld_i_fake_pseudo_local_greedy():
    check_fake('hi', pseudo=True, dowcor=True, local=True, absmean=0.400, threshold=1e-4, greedy=True)


def test_hirshfeld_i_fake_pseudo_global_greedy():
    check_fake('hi', pseudo=True, dowcor=True, local=False, absmean=0.400, threshold=1e-4, greedy=True)


def test_hirshfeld_i_fake_pseudo_nowcor_local_greedy():
    check_fake('hi', pseudo=True, dowcor=True, local=True, absmean=0.400, threshold=1e-4, greedy=True)


def test_hirshfeld_i_fake_pseudo_nowcor_global_greedy():
    check_fake('hi', pseudo=True, dowcor=True, local=False, absmean=0.400, threshold=1e-4, greedy=True)




def test_hirshfeld_e_fake_local():
    check_fake('he', pseudo=False, dowcor=True, local=True, absmean=0.323, threshold=1e-4)


def test_hirshfeld_e_fake_global():
    check_fake('he', pseudo=False, dowcor=True, local=False, absmean=0.373, threshold=1e-4)


def test_hirshfeld_e_fake_pseudo_local():
    check_fake('he', pseudo=True, dowcor=True, local=True, absmean=0.396, threshold=1e-4)


def test_hirshfeld_e_fake_pseudo_global():
    check_fake('he', pseudo=True, dowcor=True, local=False, absmean=0.396, threshold=1e-4)


def test_hirshfeld_e_fake_pseudo_nowcor_local():
    check_fake('he', pseudo=True, dowcor=True, local=True, absmean=0.396, threshold=1e-4)


def test_hirshfeld_e_fake_pseudo_nowcor_global():
    check_fake('he', pseudo=True, dowcor=True, local=False, absmean=0.396, threshold=1e-4)


def test_hirshfeld_e_fake_local_greedy():
    check_fake('he', pseudo=False, dowcor=True, local=True, absmean=0.323, threshold=1e-4, greedy=True)


def test_hirshfeld_e_fake_global_greedy():
    check_fake('he', pseudo=False, dowcor=True, local=False, absmean=0.374, threshold=1e-4, greedy=True)


def test_hirshfeld_e_fake_pseudo_local_greedy():
    check_fake('he', pseudo=True, dowcor=True, local=True, absmean=0.396, threshold=1e-4, greedy=True)


def test_hirshfeld_e_fake_pseudo_global_greedy():
    check_fake('he', pseudo=True, dowcor=True, local=False, absmean=0.396, threshold=1e-4, greedy=True)


def test_hirshfeld_e_fake_pseudo_nowcor_local_greedy():
    check_fake('he', pseudo=True, dowcor=True, local=True, absmean=0.396, threshold=1e-4, greedy=True)


def test_hirshfeld_e_fake_pseudo_nowcor_global_greedy():
    check_fake('he', pseudo=True, dowcor=True, local=False, absmean=0.396, threshold=1e-4, greedy=True)
