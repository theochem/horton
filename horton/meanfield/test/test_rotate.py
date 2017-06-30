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
from nose.plugins.attrib import attr

from horton import *  # pylint: disable=wildcard-import,unused-wildcard-import


@attr('slow')
def test_rotation_energy():
    mol = IOData.from_file(context.get_fn('test/he_spdf_orbital.fchk'))
    kin = mol.obasis.compute_kinetic()
    e0 = np.einsum('ab,ba', kin, mol.orb_alpha.to_dm())
    for irep in xrange(100):
        rmat = get_random_rotation()
        mol.orb_alpha.coeffs[:] = rotate_coeffs(mol.orb_alpha.coeffs, mol.obasis, rmat)
        e1 = np.einsum('ab,ba', kin, mol.orb_alpha.to_dm())
        assert abs(e0 - e1) < 1e-10


def test_rotation_sp():
    mol = IOData.from_file(context.get_fn('test/he_sp_orbital.fchk'))
    rmat = np.array([[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
    assert (mol.orb_alpha.coeffs[5:7,3:5] == [[0, 1], [1, 0]]).all()
    mol.orb_alpha.coeffs[:] = rotate_coeffs(mol.orb_alpha.coeffs, mol.obasis, rmat)
    assert (mol.orb_alpha.coeffs[5:7,3:5] == [[-1, 0], [0, 1]]).all()


def test_rotation_orhonormal():
    obasis = get_gobasis(np.zeros((1, 3)), np.array([10]), 'cc-pvtz', pure=False)
    overlap = obasis.compute_overlap()

    def helper(begin, end):
        # prepare random orbitals in selected range
        norb = end - begin
        assert norb > 0
        orb = Orbitals(obasis.nbasis)
        orb.occupations[:norb] = 1
        # fill with random orbitals and lowdin orthogonalize
        #orb.coeffs[begin:end,:norb] = np.random.normal(0, 1, (norb, norb))
        orb.coeffs[begin:end,:norb] = np.identity(norb)
        grammian = np.dot(orb.coeffs[:,:norb].T, np.dot(overlap, orb.coeffs[:,:norb]))
        evals, evecs = np.linalg.eigh(grammian)
        orb.coeffs[:,:norb] = np.dot(orb.coeffs[:,:norb], evecs)/evals**0.5
        assert (orb.coeffs[:begin,:norb] == 0.0).all()
        assert (orb.coeffs[end:,:norb] == 0.0).all()
        orb.check_normalization(overlap)
        # apply rotation and check normalization again
        rmat = get_random_rotation()
        #rmat = np.array([[0, 1, 0], [0, 0, 1], [1, 0, 0]])
        #rmat = get_rotation_matrix(np.array([0, 0, 1]), np.pi/4)
        orb.coeffs[:,:norb] = rotate_coeffs(orb.coeffs[:,:norb], obasis, rmat)
        orb.check_normalization(overlap)

    helper( 0,  4) # all s-type basis functions
    helper( 4, 13) # all p-type basis functions
    helper(13, 14) # just one d-type basis function
    helper(13, 19) # first set of d-type basis functions
    helper(13, 25) # all d-type basis functions
    helper(25, 35) # all f-type basis functions
