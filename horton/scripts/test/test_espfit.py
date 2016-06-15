# -*- coding: utf-8 -*-
# HORTON: Helpful Open-source Research TOol for N-fermion systems.
# Copyright (C) 2011-2016 The HORTON Development Team
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


import os, h5py as h5, numpy as np
from nose.plugins.attrib import attr

from horton import *  # pylint: disable=wildcard-import,unused-wildcard-import
from horton.test.common import check_script, tmpdir
from horton.scripts.test.common import copy_files, check_files, write_random_lta_cube
from horton.scripts.espfit import *


def test_wdens():
    assert parse_wdens('fubar.cube') == ('fubar.cube', -9, 0.8)
    assert parse_wdens('fubar.cube:-6') == ('fubar.cube', -6, 0.8)
    assert parse_wdens('fubar.cube:-6:0.5') == ('fubar.cube', -6, 0.5)

def test_wnear():
    assert parse_wnear('1:1.0') == {1: (1.0*angstrom, 0.5*angstrom)}
    assert parse_wnear('1:1.0:0.3') == {1: (1.0*angstrom, 0.3*angstrom)}
    assert parse_wnear(['1:1.0', '2:1.2']) == {1: (1.0*angstrom, 0.5*angstrom), 2: (1.2*angstrom, 0.6*angstrom)}
    assert parse_wnear(['1:1.0:0.3', '2:1.2:0.2']) == {1: (1.0*angstrom, 0.3*angstrom), 2: (1.2*angstrom, 0.2*angstrom)}

def test_wfar():
    assert parse_wfar('4.3') == (4.3*angstrom, 1.0*angstrom)
    assert parse_wfar('4.2:0.3') == (4.2*angstrom, 0.3*angstrom)


@attr('slow')
def test_scripts():
    # Generate some random system with random esp data
    natom = 5
    numbers = np.random.randint(1, 20, natom)
    coordinates = np.random.uniform(0, 10, (natom, 3))
    origin = np.zeros(3, float)
    grid_rvecs = np.identity(3, float)*1.0
    shape = np.array([10, 10, 10])
    pbc = np.ones(3, int)
    ugrid = UniformGrid(origin, grid_rvecs, shape, pbc)
    esp_cube_data = np.random.uniform(-1, 1, shape)
    rho_cube_data = np.random.uniform(-1, 1, shape)
    mol_esp = IOData(coordinates=coordinates, numbers=numbers, grid=ugrid, cube_data=esp_cube_data)
    mol_rho = mol_esp.copy()
    mol_rho.cube_data = rho_cube_data

    # Write the cube file to the tmpdir and run scripts (run 1)
    with tmpdir('horton.scripts.test.test_espfit.test_scripts') as dn:
        mol_esp.to_file(os.path.join(dn, 'esp.cube'))
        check_script('horton-esp-cost.py esp.cube esp.h5 --wnear=0:1.0:0.5', dn)
        check_script('horton-esp-fit.py esp.h5 other.h5', dn)
        check_script('horton-esp-test.py esp.h5 other.h5:charges foo.h5', dn)
        check_script('horton-esp-gen.py other.h5:charges esp.cube gen.h5', dn)
        check_files(dn, ['esp.h5', 'other.h5', 'foo.h5', 'gen.h5'])

    # Write the cube file to the tmpdir and run scripts (run 2)
    with tmpdir('horton.scripts.test.test_espfit.test_scripts2') as dn:
        mol_esp.to_file(os.path.join(dn, 'esp.cube'))
        mol_rho.to_file(os.path.join(dn, 'rho.cube'))
        check_script('horton-esp-cost.py esp.cube esp.h5 --wnear=0:1.0:0.5 --wdens=rho.cube --wsave=weight.cube', dn)
        check_files(dn, ['esp.h5', 'weight.cube'])
        check_script('horton-esp-fit.py esp.h5 other.h5', dn)
        check_script('horton-esp-test.py esp.h5 other.h5:charges foo.h5', dn)
        check_script('horton-esp-gen.py other.h5:charges esp.cube gen.h5', dn)
        check_files(dn, ['esp.h5', 'other.h5', 'foo.h5', 'gen.h5'])


@attr('slow')
def test_scripts_symmetry():
    # Write the cube file to the tmpdir and run scripts
    with tmpdir('horton.scripts.test.test_espfit.test_scripts_symmetry') as dn:
        # prepare files
        write_random_lta_cube(dn, 'esp.cube')
        copy_files(dn, ['lta_gulp.cif'])
        # run scripts
        check_script('horton-esp-cost.py esp.cube esp.h5 --wnear=0:1.0:0.5 --rcut=4 --alpha-scale=0.1', dn)
        check_files(dn, ['esp.h5'])
        check_script('horton-esp-fit.py esp.h5 other.h5 --symmetry esp.cube lta_gulp.cif', dn)
        mol_sym = IOData.from_file('%s/lta_gulp.cif' % dn)
        with h5.File(os.path.join(dn, 'other.h5')) as f:
            assert 'symmetry' in f
            assert f['symmetry/charges'].shape == (mol_sym.symmetry.natom, 2)


def test_max_at_edge():
    weights = np.array([[[0.0, 1.0], [2.0, 3.0]], [[4.0, 5.0], [6.0, 7.0]]])
    assert max_at_edge(weights, [1,1,1]) == 0.0
    assert max_at_edge(weights, [1,1,0]) == 7.0
    weights = np.array([[[0.0, 1.0, 2.0], [2.0, 3.0, 4.0]], [[4.0, 5.0, 6.0], [6.0, 9.0, 8.0]]])
    assert max_at_edge(weights, [1,1,1]) == 0.0
    assert max_at_edge(weights, [1,1,0]) == 8.0
    assert max_at_edge(weights, [1,0,1]) == 9.0
    assert max_at_edge(weights, [0,1,1]) == 9.0
