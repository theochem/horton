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


from horton import *
import numpy as np
from horton.test.common import check_delta


@attr('slow')
def test_ap1rog_one_dm():
    fn_xyz = context.get_fn('test/li2.xyz')
    mol = IOData.from_file(fn_xyz)
    obasis = get_gobasis(mol.coordinates, mol.numbers, 'cc-pvdz')
    lf = DenseLinalgFactory(obasis.nbasis)
    occ_model = AufbauOccModel(3)
    exp_alpha = lf.create_expansion(obasis.nbasis)
    olp = obasis.compute_overlap(lf)
    kin = obasis.compute_kinetic(lf)
    na = obasis.compute_nuclear_attraction(mol.coordinates, mol.pseudo_numbers, lf)
    guess_core_hamiltonian(olp, kin, na, exp_alpha)

    er = obasis.compute_electron_repulsion(lf)
    external = {'nn': compute_nucnuc(mol.coordinates, mol.pseudo_numbers)}
    terms = [
        RTwoIndexTerm(kin, 'kin'),
        RDirectTerm(er, 'hartree'),
        RExchangeTerm(er, 'x_hf'),
        RTwoIndexTerm(na, 'ne'),
    ]
    ham = REffHam(terms, external)
    scf_solver = PlainSCFSolver()
    scf_solver(ham, lf, olp, occ_model, exp_alpha)

    one = lf.create_two_index(obasis.nbasis)
    one.iadd(kin)
    one.iadd(na)

    # Do AP1roG optimization:
    geminal_solver = RAp1rog(lf, occ_model)
    energy, g, l = geminal_solver(one, er, external['nn'], exp_alpha, olp, True, **{'checkpoint': -1, 'maxiter': {'orbiter': 0}})

    one_mo_ = lf.create_two_index()
    one_mo_.assign_two_index_transform(one, exp_alpha)
    two_mo = []
    output = geminal_solver.lf.create_four_index()
    output.assign_four_index_transform(er, exp_alpha, exp_alpha, exp_alpha, exp_alpha, 'tensordot')
    two_mo.append(output)

    def fun(x):
        one_mo = []
        one_mo.append(geminal_solver.lf.create_two_index())
        one_mo[0].assign(x.reshape(28,28))

        geminal_solver.clear_auxmatrix()
        geminal_solver.update_auxmatrix('scf', two_mo, one_mo)

        iiaa = geminal_solver.get_auxmatrix('gppqq')
        iaia = geminal_solver.get_auxmatrix('lpqpq')
        fock = geminal_solver.get_auxmatrix('fock')
        one = geminal_solver.get_auxmatrix('t')
        coeff = geminal_solver.geminal._array
        lcoeff = geminal_solver.lagrange._array

        lagrangian = geminal_solver.compute_total_energy()
        lagrangian += np.dot(lcoeff.ravel(order='C'), geminal_solver.vector_function_geminal(coeff, iiaa, iaia, one, fock))
        return lagrangian


    def fun_deriv(x):
        one_mo = []
        one_mo.append(geminal_solver.lf.create_two_index())
        one_mo[0].assign(x.reshape(28,28))
        geminal_solver.clear_auxmatrix()
        geminal_solver.update_auxmatrix('scf', two_mo, one_mo)

        guesst = geminal_solver.generate_guess({'type': 'random', 'factor': -0.1})
        # Optimize OAP1roG wavefunction amplitudes:
        coeff = geminal_solver.solve_geminal(guesst, {'wfn': 'krylov'}, 10e-12, 128)

        # Optimize OAP1roG Lagrange multipliers (lambda equations):
        lcoeff = geminal_solver.solve_lagrange(guesst, {'lagrange': 'krylov'}, 10e-12, 128)
        onebody1 = geminal_solver.lf.create_two_index(3,25)
        onebody2 = geminal_solver.lf.create_two_index(3,25)
        onebody1.assign(coeff)
        onebody2.assign(lcoeff)

        onedm = geminal_solver.lf.create_one_index()
        geminal_solver.compute_1dm(onedm, onebody1, onebody2, factor=2.0)
        a = np.zeros((28,28))
        np.fill_diagonal(a, onedm._array.T)
        return a.ravel()

    x = one_mo_._array.ravel()
    dxs = np.random.rand(100, 28*28)*0.00001
    check_delta(fun, fun_deriv, x, dxs)
