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

def test_gradient_ap1rog_cs_scf():
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
    energy, g, l = geminal_solver(one, er, external['nn'], exp_alpha, olp, True, **{'solver': {'wfn': 'hybr', 'lagrange': 'hybr'}, 'checkpoint': -1, 'maxiter': {'orbiter': 20}})

    x = np.zeros(378)
    dxs = np.random.rand(500, 378)*0.001
    check_delta(x, dxs, geminal_solver, exp_alpha, one, er)

def fun(x, ham, orb, one, two):
    tmp = orb.copy()
    rotation = ham.get_rotation_matrix(x)
    ham.rotate_orbitals(tmp, rotation)
    ham.clear_matrix()
    one_mo, two_mo = ham.update_mo_integrals(one, two, 'tensordot', tmp)

    ham.update_matrix('all-scf', two_mo, one_mo)

    kmatrix = ham.get_matrix('vpp')
    vmatrix = ham.get_matrix('v')
    vmatrixa = ham.get_matrix('va')
    one_mo = ham.get_matrix('t')
    energy = ham.get_reference_energy()
    coeff = ham.geminal._array
    lcoeff = ham.lagrange._array
    lagrangian = ham.get_total_energy()
    lagrangian += np.dot(lcoeff.ravel(order='C'), ham.vector_function_ap1rog(coeff,kmatrix,vmatrix,vmatrixa,one_mo,energy))
    return lagrangian

def fun_deriv(x, ham, orb, one, two):
    tmp = orb.copy()
    rotation = ham.get_rotation_matrix(x)
    ham.rotate_orbitals(tmp, rotation)
    ham.clear_matrix()
    one_mo, two_mo = ham.update_mo_integrals(one, two, 'tensordot', tmp)

    ham.update_matrix('all-scf', two_mo, one_mo)

    gradient = ham.get_orbital_gradient()
    return gradient

def check_delta(x, dxs, ham, orb, one, two):
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
    f0 = fun(x, ham, orb, one, two)
    grad0 = fun_deriv(x, ham, orb, one, two)
    for dx in dxs:
        rot = ham.get_rotation_matrix(dx)
        tril = np.tril_indices(28,-1)
        f1 = fun(x+dx, ham, orb, one, two)
        grad1 = fun_deriv(x+dx, ham, orb, one, two)
        grad = 0.5*(grad0+grad1)
        d1 = f1 - f0
        if hasattr(d1, '__iter__'):
            norm = np.linalg.norm
        else:
            norm = abs
        d2 = np.dot(grad.ravel(), rot[tril].ravel())

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
