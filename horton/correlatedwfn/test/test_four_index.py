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


import numpy as np
from nose.tools import assert_raises
from horton import *


def test_tensordot_rs():
    fn_xyz = context.get_fn('test/h2.xyz')
    mol = Molecule.from_file(fn_xyz)
    obasis = get_gobasis(mol.coordinates, mol.numbers, 'cc-pvdz')
    lf = DenseLinalgFactory(obasis.nbasis)
    occ_model = AufbauOccModel(1)
    exp_alpha = lf.create_expansion(obasis.nbasis)
    olp = obasis.compute_overlap(lf)
    kin = obasis.compute_kinetic(lf)
    na = obasis.compute_nuclear_attraction(mol.coordinates, mol.pseudo_numbers, lf)
    guess_core_hamiltonian(olp, kin, na, exp_alpha)

    er = obasis.compute_electron_repulsion(lf)
    external = {'nn': compute_nucnuc(mol.coordinates, mol.pseudo_numbers)}
    terms = [
        ROneBodyTerm(kin, 'kin'),
        RDirectTerm(er, 'hartree'),
        RExchangeTerm(er, 'x_hf'),
        ROneBodyTerm(na, 'ne'),
    ]
    ham = REffHam(terms, external)
    scf_solver = PlainSCFSolver()
    scf_solver(ham, lf, olp, occ_model, exp_alpha)

    one = lf.create_one_body(obasis.nbasis)
    one.iadd(kin)
    one.iadd(na)

    # Do AP1roG optimization:
    geminal_solver = RAp1rog(lf, occ_model)

    one_mo, two_mo = geminal_solver.update_mo_integrals(one, er, 'tensordot', exp_alpha)

    assert max(abs((one_mo[0]._array.ravel())-(one_mo[0]._array.transpose(1,0).ravel()))) < 1e-12

    assert max(abs((two_mo[0]._array.ravel())-(two_mo[0]._array.transpose(1,0,3,2).ravel()))) < 1e-12
    assert max(abs((two_mo[0]._array.ravel())-(two_mo[0]._array.transpose(2,3,0,1).ravel()))) < 1e-12
    assert max(abs((two_mo[0]._array.ravel())-(two_mo[0]._array.transpose(3,2,1,0).ravel()))) < 1e-12
    assert max(abs((two_mo[0]._array.ravel())-(two_mo[0]._array.transpose(2,1,0,3).ravel()))) < 1e-12
    assert max(abs((two_mo[0]._array.ravel())-(two_mo[0]._array.transpose(3,0,1,2).ravel()))) < 1e-12
    assert max(abs((two_mo[0]._array.ravel())-(two_mo[0]._array.transpose(0,3,2,1).ravel()))) < 1e-12
    assert max(abs((two_mo[0]._array.ravel())-(two_mo[0]._array.transpose(1,2,3,0).ravel()))) < 1e-12

def test_einstein_rs():
    fn_xyz = context.get_fn('test/h2.xyz')
    mol = Molecule.from_file(fn_xyz)
    obasis = get_gobasis(mol.coordinates, mol.numbers, 'cc-pvdz')
    lf = DenseLinalgFactory(obasis.nbasis)
    occ_model = AufbauOccModel(1)
    exp_alpha = lf.create_expansion(obasis.nbasis)
    olp = obasis.compute_overlap(lf)
    kin = obasis.compute_kinetic(lf)
    na = obasis.compute_nuclear_attraction(mol.coordinates, mol.pseudo_numbers, lf)
    guess_core_hamiltonian(olp, kin, na, exp_alpha)

    er = obasis.compute_electron_repulsion(lf)
    external = {'nn': compute_nucnuc(mol.coordinates, mol.pseudo_numbers)}
    terms = [
        ROneBodyTerm(kin, 'kin'),
        RDirectTerm(er, 'hartree'),
        RExchangeTerm(er, 'x_hf'),
        ROneBodyTerm(na, 'ne'),
    ]
    ham = REffHam(terms, external)
    scf_solver = PlainSCFSolver()
    scf_solver(ham, lf, olp, occ_model, exp_alpha)

    one = lf.create_one_body(obasis.nbasis)
    one.iadd(kin)
    one.iadd(na)

    # Do AP1roG optimization:
    geminal_solver = RAp1rog(lf, occ_model)

    one_mo, two_mo = geminal_solver.update_mo_integrals(one, er, 'einstein', exp_alpha)

    assert max(abs((one_mo[0]._array.ravel())-(one_mo[0]._array.transpose(1,0).ravel()))) < 1e-12

    assert max(abs((two_mo[0]._array.ravel())-(two_mo[0]._array.transpose(1,0,3,2).ravel()))) < 1e-12
    assert max(abs((two_mo[0]._array.ravel())-(two_mo[0]._array.transpose(2,3,0,1).ravel()))) < 1e-12
    assert max(abs((two_mo[0]._array.ravel())-(two_mo[0]._array.transpose(3,2,1,0).ravel()))) < 1e-12
    assert max(abs((two_mo[0]._array.ravel())-(two_mo[0]._array.transpose(2,1,0,3).ravel()))) < 1e-12
    assert max(abs((two_mo[0]._array.ravel())-(two_mo[0]._array.transpose(3,0,1,2).ravel()))) < 1e-12
    assert max(abs((two_mo[0]._array.ravel())-(two_mo[0]._array.transpose(0,3,2,1).ravel()))) < 1e-12
    assert max(abs((two_mo[0]._array.ravel())-(two_mo[0]._array.transpose(1,2,3,0).ravel()))) < 1e-12
