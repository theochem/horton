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


def test_rotation_energy():
    mol = Molecule.from_file(context.get_fn('test/he_spdf_orbital.fchk'))
    kin = mol.obasis.compute_kinetic(mol.lf)
    e0 = kin.expectation_value(mol.exp_alpha.to_dm())
    for irep in xrange(100):
        rmat = get_random_rotation()
        mol.exp_alpha.coeffs[:] = rotate_coeffs(mol.exp_alpha.coeffs, mol.obasis, rmat)
        e1 = kin.expectation_value(mol.exp_alpha.to_dm())
        assert abs(e0 - e1) < 1e-10


def test_rotation_sp():
    mol = Molecule.from_file(context.get_fn('test/he_sp_orbital.fchk'))
    rmat = np.array([[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
    assert (mol.exp_alpha.coeffs[5:7,3:5] == [[0, 1], [1, 0]]).all()
    mol.exp_alpha.coeffs[:] = rotate_coeffs(mol.exp_alpha.coeffs, mol.obasis, rmat)
    assert (mol.exp_alpha.coeffs[5:7,3:5] == [[-1, 0], [0, 1]]).all()
