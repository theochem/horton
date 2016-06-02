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


import numpy as np
from nose.plugins.attrib import attr

from horton import *
from horton.test.common import tmpdir


def prepare_mol(basis, scale=1.0):
    fn_xyz = context.get_fn('test/h2.xyz')
    mol = IOData.from_file(fn_xyz)
    mol.coordinates *= scale
    obasis = get_gobasis(mol.coordinates, mol.numbers, basis)
    lf = DenseLinalgFactory(obasis.nbasis)
    olp = obasis.compute_overlap(lf)
    kin = obasis.compute_kinetic(lf)
    na = obasis.compute_nuclear_attraction(mol.coordinates, mol.pseudo_numbers, lf)
    one = kin.copy()
    one.iadd(na)
    er = obasis.compute_electron_repulsion(lf)

    return mol, lf, olp, one, er


def prepare_hf(basis, scale=1.0):
    mol, lf, olp, one, er = prepare_mol(basis, scale)

    exp_alpha = lf.create_expansion()
    guess_core_hamiltonian(olp, one, exp_alpha)
    external = {'nn': compute_nucnuc(mol.coordinates, mol.pseudo_numbers)}
    terms = [
        RTwoIndexTerm(one, 'one'),
        RDirectTerm(er, 'hartree'),
        RExchangeTerm(er, 'x_hf'),
    ]
    ham = REffHam(terms, external)
    occ_model = AufbauOccModel(1)
    scf_solver = PlainSCFSolver()
    scf_solver(ham, lf, olp, occ_model, exp_alpha)

    return lf, occ_model, one, er, external, exp_alpha, olp


def test_ap1rog_cs():
    lf, occ_model, one, er, external, exp_alpha, olp = prepare_hf('6-31G')

    # Do AP1roG optimization:
    geminal_solver = RAp1rog(lf, occ_model)
    guess = np.array([-0.08, -0.05, -0.03])
    energy, g = geminal_solver(one, er, external['nn'], exp_alpha, olp, False, **{'guess': {'geminal': guess}})
    assert (abs(energy - -1.143420629378) < 1e-6)


@attr('slow')
def test_ap1rog_cs_scf():
    lf, occ_model, one, er, external, exp_alpha, olp = prepare_hf('6-31G')

    # Do AP1roG optimization:
    geminal_solver = RAp1rog(lf, occ_model)
    guess = np.array([-0.1, -0.05, -0.02])
    energy, g, l = geminal_solver(one, er, external['nn'], exp_alpha, olp, True, **{'checkpoint': -1, 'guess': {'geminal': guess}})
    assert (abs(energy - -1.151686291339) < 1e-6)


@attr('slow')
def test_ap1rog_cs_scf_restart():
    lf, occ_model, one, er, external, exp_alpha, olp = prepare_hf('6-31G')

    # Do AP1roG optimization:
    geminal_solver = RAp1rog(lf, occ_model)
    guess = np.array([-0.11, -0.05, -0.02])

    with tmpdir('horton.correlatedwfn.test.test_geminals.test_ap1rog_cs_scf_restart') as dn:
        checkpoint_fn = '%s/checkpoint.h5' % dn
        energy, g, l = geminal_solver(one, er, external['nn'], exp_alpha, olp,
                                      True, guess={'geminal': guess},
                                      checkpoint_fn=checkpoint_fn)
        assert (abs(energy - -1.151686291339) < 1e-6)

        old = IOData.from_file(checkpoint_fn)

    assert hasattr(old, 'olp')
    assert hasattr(old, 'exp_alpha')

    # Update to slightly stretch geometry of H2
    mol1, lf, olp1, one1, er1 = prepare_mol('6-31G', 1.05)

    # re-orthogonalize orbitals
    project_orbitals_ortho(old.olp, olp1, old.exp_alpha, exp_alpha)

    # recompute
    energy, c, l = geminal_solver(one1, er1, external['nn'], exp_alpha, olp1,
                                  True, guess={'geminal': guess},
                                  checkpoint=-1)

    assert (abs(energy - -1.11221603918) < 1e-6)
