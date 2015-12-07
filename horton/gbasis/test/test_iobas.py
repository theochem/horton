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


import numpy as np

from horton import *


def test_str_to_shell_types_cart():
    assert str_to_shell_types('s') == [0]
    assert str_to_shell_types('S') == [0]
    assert str_to_shell_types('Ss') == [0,0]
    assert str_to_shell_types('SP') == [0,1]
    assert str_to_shell_types('SDD') == [0,2,2]


def test_str_to_shell_types_pure():
    assert str_to_shell_types('s', True) == [0]
    assert str_to_shell_types('S', True) == [0]
    assert str_to_shell_types('Ss', True) == [0,0]
    assert str_to_shell_types('SP', True) == [0,1]
    assert str_to_shell_types('SDF', True) == [0,-2,-3]


def test_shell_type_to_str():
    assert shell_type_to_str(0) == 's'
    assert shell_type_to_str(-1) == 'p'
    assert shell_type_to_str(1) == 'p'
    assert shell_type_to_str(-3) == 'f'
    assert shell_type_to_str(4) == 'g'


def test_fortran_float():
    assert fortran_float('1e-2') == fortran_float('1d-2')
    assert float('1e-2') == fortran_float('1d-2')
    assert fortran_float('1E-2') == fortran_float('1D-2')
    assert float('1E-2') == fortran_float('1D-2')


def test_go_basis_desc_neon_sto3g():
    gobasis_set = []
    gobasis_set.append(GOBasisFamily('STO-3G', filename=context.get_fn('basis/sto-3g.nwchem')))
    gobasis_set.append(GOBasisFamily('STO-3G', filename=context.get_fn('basis/sto-3g.gbs')))
    for gobasis in gobasis_set:
        obasis = get_gobasis(np.array([[0.0,0.0,0.0]]), np.array([2]), gobasis)
        assert (obasis.shell_map == np.array([0])).all()
        assert (obasis.nprims == np.array([3])).all()
        assert (obasis.shell_types == np.array([0])).all()
        assert abs(obasis.alphas - np.array([6.36242139, 1.15892300, 0.31364979])).all() < 1e-10
        assert abs(obasis.con_coeffs - np.array([0.15432897, 0.53532814, 0.44463454])).all() < 1e-10


def test_go_basis_desc_hydrogen_321g():
    gobasis_set = []
    gobasis_set.append(GOBasisFamily('3-21G', filename=context.get_fn('basis/3-21g.nwchem')))
    gobasis_set.append(GOBasisFamily('3-21G', filename=context.get_fn('basis/3-21g.gbs')))
    for gobasis in gobasis_set:
        obasis = get_gobasis(np.array([[0.0,0.0,0.0]]), np.array([1]), gobasis)
        assert (obasis.shell_map == np.array([0,0])).all()
        assert (obasis.nprims == np.array([2,1])).all()
        assert (obasis.shell_types == np.array([0,0])).all()
        assert abs(obasis.alphas - np.array([5.4471780, 0.8245470, 0.1831920])).all() < 1e-10
        assert abs(obasis.con_coeffs - np.array([0.1562850, 0.9046910, 1.0000000])).all() < 1e-10
        assert obasis.nbasis == 2


def test_go_basis_family_lithium_321g():
    bf = go_basis_families['3-21g']
    assert len(bf.get(3).bcs) == 5


def test_go_basis_desc_lithium_321g():
    gobasis_set = []
    gobasis_set.append(GOBasisFamily('3-21G', filename=context.get_fn('basis/3-21g.nwchem')))
    gobasis_set.append(GOBasisFamily('3-21G', filename=context.get_fn('basis/3-21g.gbs')))
    for gobasis in gobasis_set:
        obasis = get_gobasis(np.array([[0.0,0.0,0.0]]), np.array([3]), gobasis)
        assert (obasis.shell_map == np.array([0,0,0,0,0])).all()
        assert (obasis.nprims == np.array([3,2,2,1,1])).all()
        assert (obasis.shell_types == np.array([0,0,1,0,1])).all()
        assert abs(obasis.alphas - np.array([
            36.8382000, 5.4817200, 1.1132700,
            0.5402050, 0.1022550, 0.5402050, 0.1022550,
            0.0285650, 0.0285650,
        ])).all() < 1e-10
        assert abs(obasis.con_coeffs - np.array([
            0.0696686, 0.3813460, 0.6817020,
            -0.2631270, 0.1615460, 1.1433900, 0.9156630,
            1.0000000, 1.0000000
        ])).all() < 1e-10
        assert obasis.nbasis == 9


def test_go_basis_desc_water_sto3g():
    fn = context.get_fn('test/water_element.xyz')
    mol = IOData.from_file(fn)
    gobasis_set = []
    gobasis_set.append(GOBasisFamily('STO-3G', filename=context.get_fn('basis/sto-3g.nwchem')))
    gobasis_set.append(GOBasisFamily('STO-3G', filename=context.get_fn('basis/sto-3g.gbs')))
    for gobasis in gobasis_set:
        obasis = get_gobasis(mol.coordinates, mol.numbers, gobasis)
        assert (obasis.shell_map == np.array([0,1,1,1,2])).all()
        assert (obasis.nprims == np.array([3,3,3,3,3])).all()
        assert (obasis.shell_types == np.array([0,0,0,1,0])).all()
        assert abs(obasis.alphas - np.array([
            3.42525091, 0.62391373, 0.16885540,
            130.7093200, 23.8088610, 6.4436083,
            5.0331513, 1.1695961, 0.3803890,
            5.0331513, 1.1695961, 0.3803890,
            3.42525091, 0.62391373, 0.16885540,
        ])).all() < 1e-10
        assert abs(obasis.con_coeffs - np.array([
            0.15432897, 0.53532814, 0.44463454,
            0.15432897, 0.53532814, 0.44463454,
            -0.09996723, 0.15591627, 0.39951283,
            0.60768372, 0.70011547, 0.39195739,
            0.15432897, 0.53532814, 0.44463454,
        ])).all() < 1e-10
        assert obasis.nbasis == 7
