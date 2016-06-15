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
# --


import numpy as np
from nose.tools import assert_raises
from nose.plugins.attrib import attr

from horton import *  # pylint: disable=wildcard-import,unused-wildcard-import

from horton.test.common import numpy_seed


@attr('slow')
def test_ap1rog_lagrange():
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
    with numpy_seed():
        energy, g = geminal_solver(one, er, external['nn'], exp_alpha, olp, False)
        geminal_solver.lagrange.assign(np.random.rand(3,25))
        x = geminal_solver.geminal._array.ravel(order='C')
        dxs = np.random.rand(200, 3*25)*(0.001)
        check_delta(x, dxs, geminal_solver)


def fun(x, ham):
    iiaa = ham.get_auxmatrix('gppqq')
    iaia = ham.get_auxmatrix('lpqpq')
    fock = ham.get_auxmatrix('fock')
    one = ham.get_auxmatrix('t')
    coeff = ham.lagrange._array

    lagrangian = ham.compute_total_energy(x.reshape(3,25))
    lagrangian += np.dot(coeff.ravel(order='C'), ham.vector_function_geminal(x.ravel(order='C'), iiaa, iaia, one, fock))
    return lagrangian


def fun_deriv(x, ham):
    iiaa = ham.get_auxmatrix('gppqq')
    iaia = ham.get_auxmatrix('lpqpq')
    fock = ham.get_auxmatrix('fock')
    one = ham.get_auxmatrix('t')
    coeff = ham.lagrange._array.ravel(order='C')
    gmat = ham.lf.create_two_index(3,25)
    gmat.assign(x)

    gradient = ham.vector_function_lagrange(coeff,gmat, iiaa, iaia, one, fock)
    return gradient.ravel(order='C')


def check_delta(x, dxs, ham):
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
    dn1s = []
    dn2s = []
    dnds = []
    f0 = fun(x, ham)
    grad0 = fun_deriv(x, ham)
    for dx in dxs:
        f1 = fun(x+dx, ham)
        grad1 = fun_deriv(x+dx, ham)
        grad = 0.5*(grad0+grad1)
        d1 = f1 - f0
        if hasattr(d1, '__iter__'):
            norm = np.linalg.norm
        else:
            norm = abs
        d2 = np.dot(grad.ravel(), dx.ravel())

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
