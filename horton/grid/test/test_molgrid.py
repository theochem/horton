# -*- coding: utf-8 -*-
# Horton is a Density Functional Theory program.
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

from horton import *



def test_integrate_hydrogen_single_1s():
    numbers = np.array([1], int)
    coordinates = np.array([[0.0, 0.0, -0.5]], float)
    sys = System(coordinates, numbers)
    int1d = TrapezoidIntegrator1D()
    rtf = ExpRTransform(1e-3, 1e1, 100)
    atspecs = (rtf, int1d, 110)

    mg = BeckeMolGrid(sys, atspecs, random_rotate=False)
    dist0 = np.sqrt(((coordinates[0] - mg.points)**2).sum(axis=1))
    fn = np.exp(-2*dist0)/np.pi
    occupation = mg.integrate(fn)
    assert abs(occupation - 1.0) < 1e-3


def test_integrate_hydrogen_pair_1s():
    numbers = np.array([1, 1], int)
    coordinates = np.array([[0.0, 0.0, -0.5], [0.0, 0.0, 0.5]], float)
    sys = System(coordinates, numbers)
    int1d = TrapezoidIntegrator1D()
    rtf = ExpRTransform(1e-3, 1e1, 100)
    atspecs = (rtf, int1d, 110)

    mg = BeckeMolGrid(sys, atspecs, random_rotate=False)
    dist0 = np.sqrt(((coordinates[0] - mg.points)**2).sum(axis=1))
    dist1 = np.sqrt(((coordinates[1] - mg.points)**2).sum(axis=1))
    fn = np.exp(-2*dist0)/np.pi + np.exp(-2*dist1)/np.pi
    occupation = mg.integrate(fn)
    assert abs(occupation - 2.0) < 1e-3


def test_integrate_hydrogen_trimer_1s():
    numbers = np.array([1, 1, 1], int)
    coordinates = np.array([[0.0, 0.0, -0.5], [0.0, 0.0, 0.5], [0.0, 0.5, 0.0]], float)
    sys = System(coordinates, numbers)
    int1d = TrapezoidIntegrator1D()
    rtf = ExpRTransform(1e-3, 1e1, 100)
    atspecs = (rtf, int1d, 110)

    mg = BeckeMolGrid(sys, atspecs, random_rotate=False)
    dist0 = np.sqrt(((coordinates[0] - mg.points)**2).sum(axis=1))
    dist1 = np.sqrt(((coordinates[1] - mg.points)**2).sum(axis=1))
    dist2 = np.sqrt(((coordinates[2] - mg.points)**2).sum(axis=1))
    fn = np.exp(-2*dist0)/np.pi + np.exp(-2*dist1)/np.pi + np.exp(-2*dist2)/np.pi
    occupation = mg.integrate(fn)
    assert abs(occupation - 3.0) < 1e-3


def test_molgrid_attrs_subgrid():
    numbers = np.array([6, 8], int)
    coordinates = np.array([[0.0, 0.2, -0.5], [0.1, 0.0, 0.5]], float)
    sys = System(coordinates, numbers)
    int1d = TrapezoidIntegrator1D()
    rtf = ExpRTransform(1e-3, 1e1, 100)
    atspecs = (rtf, int1d, 110)
    mg = BeckeMolGrid(sys, atspecs, keep_subgrids=True)

    assert mg.size == 2*110*100
    assert mg.points.shape == (mg.size, 3)
    assert mg.weights.shape == (mg.size,)
    assert len(mg.subgrids) == 2
    assert mg.system == sys
    assert len(mg.atspecs) == 2
    assert mg.k == 3
    assert mg.random_rotate

    for i in xrange(2):
        atspec = mg.atspecs[i]
        assert atspec[0] == rtf
        assert atspec[1] == int1d
        assert (atspec[2] == np.array([110]*100)).all()

        atgrid = mg.subgrids[i]
        assert isinstance(atgrid, AtomicGrid)
        assert atgrid.size == 100*110
        assert atgrid.points.shape == (100*110, 3)
        assert atgrid.weights.shape == (100*110,)
        assert atgrid.av_weights.shape == (100*110,)
        assert atgrid.subgrids is None
        assert atgrid.number == numbers[i]
        assert (atgrid.center == coordinates[i]).all()
        assert atgrid.rtransform == rtf
        assert (atgrid.nlls == [110]*100).all()
        assert atgrid.nsphere == 100
        assert atgrid.random_rotate


def test_molgrid_attrs():
    numbers = np.array([6, 8], int)
    coordinates = np.array([[0.0, 0.2, -0.5], [0.1, 0.0, 0.5]], float)
    sys = System(coordinates, numbers)
    int1d = TrapezoidIntegrator1D()
    rtf = ExpRTransform(1e-3, 1e1, 100)
    atspecs = (rtf, int1d, 110)
    mg = BeckeMolGrid(sys, atspecs, keep_subgrids=False)

    assert mg.size == 2*110*100
    assert mg.points.shape == (mg.size, 3)
    assert mg.weights.shape == (mg.size,)
    assert mg.subgrids is None
    assert mg.system == sys
    assert len(mg.atspecs) == 2
    assert mg.k == 3
    assert mg.random_rotate


def test_custom_grid_term():
    fn_fchk = context.get_fn('test/n2_hfs_sto3g.fchk')
    sys = System.from_file(fn_fchk)
    int1d = SimpsonIntegrator1D()
    rtf = ExpRTransform(1e-3, 1e1, 100)
    grid = BeckeMolGrid(sys, (rtf, int1d, 110), random_rotate=False)

    # Without perturbation
    ham = Hamiltonian(sys, [HartreeFock()])
    assert convergence_error(ham) > 1e-8
    assert converge_scf(ham)
    assert convergence_error(ham) < 1e-8
    energy0 = ham.compute_energy()

    # Construct some becke weights for the first atom and use it as a potential.
    potential = np.ones(grid.size, float)
    radii = np.ones(sys.natom, float)
    becke_helper_atom(grid.points, potential, radii, sys.coordinates, 0, 3)

    # Apply the perturbation with oposite signs and check that, because of
    # symmetry, the energy of the perturbed wavefunction is the same in both
    # cases, and higher than the unperturbed.
    energy1_old = None
    for scale in 0.1, -0.1:
        # With perturbation
        perturbation = CustomGridFixedTerm(grid, scale*potential, 'pert')
        ham = Hamiltonian(sys, [HartreeFock(), perturbation])
        assert convergence_error(ham) > 1e-8
        assert converge_scf_oda(ham)
        assert convergence_error(ham) < 1e-8
        energy1 = ham.compute_energy()
        energy1 -= sys.props['energy_pert']

        assert energy1 > energy0
        if energy1_old is None:
            energy1_old = energy1
        else:
            assert abs(energy1 - energy1_old) < 1e-7


def test_family():
    numbers = np.array([6, 8], int)
    coordinates = np.array([[0.0, 0.2, -0.5], [0.1, 0.0, 0.5]], float)
    sys = System(coordinates, numbers)
    grid = BeckeMolGrid(sys, random_rotate=False)
    assert grid.size == 1434+1300
