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

def test_ap1rog_two_dm():
    fn_xyz = context.get_fn('test/li2.xyz')
    mol = Molecule.from_file(fn_xyz)
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
    energy, g, l = geminal_solver(one, er, external['nn'], exp_alpha, olp, True, **{'solver': {'wfn': 'hybr', 'lagrange': 'hybr'}, 'checkpoint': -1, 'maxiter': {'orbiter': 10}})

    one_mo = []
    one_mo.append(lf.create_two_index())
    one_mo[0].apply_2index_trans(one, exp_alpha)
    two_mo_ = []
    two_mo_.append(geminal_solver.apply_4index_trans(er, exp_alpha, exp_alpha, exp_alpha, exp_alpha))

    x = two_mo_[0]._array
    dxs = np.random.rand(30, 28, 28, 28, 28)*0.000001
    check_delta(x, dxs, geminal_solver, one_mo)

def fun(x, ham, one_mo):
    two_mo = []
    two_mo.append(ham.lf.create_four_index())
    two_mo[0].assign_array(x)

    ham.clear_matrix()
    ham.update_matrix('all-scf', two_mo, one_mo)

    kmatrix = ham.get_matrix('vpp')
    vmatrix = ham.get_matrix('v')
    vmatrixa = ham.get_matrix('va')
    one = ham.get_matrix('t')
    energy = ham.get_reference_energy()
    coeff = ham.geminal._array
    lcoeff = ham.lagrange._array
    lagrangian = ham.get_total_energy()
    lagrangian += np.dot(lcoeff.ravel(order='C'), ham.vector_function_ap1rog(coeff,kmatrix,vmatrix,vmatrixa,one,energy))
    return lagrangian

def fun_deriv(x, ham, one_mo):
    two_mo = []
    two_mo.append(ham.lf.create_four_index())
    two_mo[0].assign_array(x)
    ham.clear_matrix()
    ham.update_matrix('all-scf', two_mo, one_mo)

    guesst = ham.generate_guess(guess={'guess': 'random', 'factor': -0.1})
    # Optimize OAP1roG wavefunction amplitudes:
    coeff = ham.solve_for_wavefunction(guesst[0], {'wfn': 'hybr'}, 10e-12, 128)

    # Optimize OAP1roG Lagrange multipliers (lambda equations):
    lcoeff = ham.solve_for_lagrange(guesst[0], {'lagrange': 'hybr'}, 10e-12, 128)
    twoindex1 = ham.lf.create_two_index(3,25)
    twoindex2 = ham.lf.create_two_index(3,25)
    twoindex1.assign_array(coeff, 3, 25)
    twoindex2.assign_array(lcoeff, 3,25)

    onedm = ham.lf.create_one_index()
    onedm.compute_response_one_dm_ap1rog(twoindex1, twoindex2, factor=1.0)
    twodm1 = ham.lf.create_two_index()
    twodm1.compute_response_two_dm_ap1rog(onedm, twoindex1, twoindex2, 'ppqq')
    twodm2 = ham.lf.create_two_index()
    twodm2.compute_response_two_dm_ap1rog(onedm, twoindex1, twoindex2, 'pqpq')
    a = np.zeros((28,28,28,28))
    for i in range(28):
        for j in range(28):
            a[i,i,j,j] += twodm1._array[i,j]
            a[i,j,i,j] += twodm2._array[i,j]
            a[i,j,i,j] += twodm2._array[i,j]
            a[i,j,j,i] -= twodm2._array[i,j]
    return a

def check_delta(x, dxs, ham, orb):
    """Check the difference between two function values using the analytical gradient

       Arguments:

       fun
            The function whose derivatives must be to be tested

       fun_deriv
            The implementation of the analytical derivatives

       x
            The argument for the reference point.

       dxs
            A list with small relative changes to x

       For every displacement in ``dxs``, the following computation is repeated:

       1) D1 = 'fun(x+dx) - fun(x)' is computed.
       2) D2 = '0.5 (fun_deriv(x+dx) + fun_deriv(x)) . dx' is computed.

       A threshold is set to the median of the D1 set. For each case where |D1|
       is larger than the threshold, |D1 - D2|, should be smaller than the
       threshold.
       """
#   if len(dxs) < 20:
#       raise ValueError('At least 20 displacements are needed for good statistics.')

    dn1s = []
    dn2s = []
    dnds = []
    f0 = fun(x, ham, orb)
    grad0 = fun_deriv(x, ham, orb)
    for dx in dxs:
        f1 = fun(x+dx, ham, orb)
        grad1 = fun_deriv(x+dx, ham, orb)
        grad = 0.5*(grad0+grad1)
        d1 = f1 - f0
        if hasattr(d1, '__iter__'):
            norm = np.linalg.norm
        else:
            norm = abs
        d2 = np.dot(grad.ravel(), dx.ravel())
        print d1, d2

        dn1s.append(norm(d1))
        dn2s.append(norm(d2))
        dnds.append(norm(d1-d2))
    dn1s = np.array(dn1s)
    dn2s = np.array(dn2s)
    dnds = np.array(dnds)

    # Get the threshold (and mask)
    threshold = np.median(dn1s)
    mask = dn1s > threshold
    # Make sure that all cases for which dn1 is above the treshold, dnd is below
    # the threshold
    if not (dnds[mask] < threshold).all():
        raise AssertionError((
            'The first order approximation on the difference is too wrong. The '
            'threshold is %.1e.\n\nDifferences:\n%s\n\nFirst order '
            'approximation to differences:\n%s\n\nAbsolute errors:\n%s')
            % (threshold,
            ' '.join('%.1e' % v for v in dn1s[mask]),
            ' '.join('%.1e' % v for v in dn2s[mask]),
            ' '.join('%.1e' % v for v in dnds[mask])
        ))
