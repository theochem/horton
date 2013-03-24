#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Horton is a Density Functional Theory program.
# Copyright (C) 2011-2012 Toon Verstraelen <Toon.Verstraelen@UGent.be>
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


import sys, argparse, os

import h5py as h5, numpy as np
from horton import System, UniformIntGrid, log, Cell, angstrom, compute_esp_grid_cube
from horton.scripts.common import parse_h5


def parse_args():
    parser = argparse.ArgumentParser(prog='horton-esp-gen.py',
        description='Generate electrostatic potential grid data from charges.')

    parser.add_argument('h5_charges',
        help='[HDF5 filename]:[HDF5 group] The HDF5 file with the charges. '
             'The results of the test will be written in a subgroup "espgrid". '
             'The HDF5 file must also contain a system group with coordinates '
             'and a cell group.')
    parser.add_argument('--qtot', '-q', default=None, type=float,
        help='The total charge of the system. When given, the charges from the '
             'HDF5 file are corrected.')

    parser.add_argument('--spacing', '-s', default=0.1, type=float,
        help='The target grid spacing in angstrom.')

    parser.add_argument('--rcut', default=10.0, type=float,
        help='The real-space cutoff for the electrostatic interactions in '
             'angstrom.')
    parser.add_argument('--alpha-scale', default=3.0, type=float,
        help='The alpha scale (alpha = alpha_scale/rcut) for the separation '
             'between short-range and long-range electrostatic interactions.')
    parser.add_argument('--gcut-scale', default=1.1, type=float,
        help='The gcut scale (gcut = gcut_scale*alpha) for the reciprocal '
             'space constribution to the electrostatic interactions.')

    return parser.parse_args()


def main():
    args = parse_args()

    # Load the charges from the HDF5 file
    fn_h5, grp_name = parse_h5(args.h5_charges)
    with h5.File(fn_h5, 'r') as f:
        grp = f[grp_name]
        charges = grp['charges'][:]
        cell = Cell.from_hdf5(f['system/cell'], None)
        sys = System(f['system/coordinates'][:], f['system/numbers'][:], cell=cell)

    if sys.cell.nvec != 3:
        raise NotImplementedError('Only 3D periodic cells are suppported.')

    # Fix total charge if requested
    if args.qtot is not None:
        charges -= (charges.sum() - args.qtot)/len(charges)

    # Store parameters in output
    results = {}
    results['qtot'] = charges.sum()
    results['spacing'] = args.spacing*angstrom
    results['rcut'] = args.rcut*angstrom
    results['alpha_scale'] = args.alpha_scale
    results['gcut_scale'] = args.gcut_scale

    # Determine the grid specification
    grid_rvecs = cell.rvecs.copy()
    shape = np.zeros(3, int)
    lengths, angles = cell.parameters
    for i in xrange(3):
        shape[i] = int(np.round(lengths[i]/(args.spacing*angstrom)))
        grid_rvecs[i] /= shape[i]
    grid_cell = Cell(grid_rvecs)
    origin = np.zeros(3, float)
    pbc = np.ones(3, int)
    ui_grid = UniformIntGrid(origin, grid_rvecs, shape, pbc)

    # Ewald parameters
    rcut = args.rcut*angstrom      # <- TODO: three lines into one routine
    alpha = args.alpha_scale/rcut  # <-
    gcut = args.gcut_scale*alpha   # <-

    # Some screen info
    if log.do_medium:
        log('Important parameters:')
        log.hline()
        log('Number of grid points:   %12i' % np.product(shape))
        log('Grid shape:                 [%8i, %8i, %8i]' % tuple(shape))
        log.hline()
        # TODO: add ewald parameters and summation ranges
        log('Computing ESP (may take a while)')

    # Allocate and compute ESP grid
    esp = np.zeros(shape, float)
    compute_esp_grid_cube(ui_grid, esp, sys.coordinates, charges, rcut, alpha, gcut)
    results['esp'] = esp

    # Store the results in an HDF5 file
    with h5.File(fn_h5) as f:
        # Get the group for the output
        grp = f[grp_name]
        if 'espgrid' in grp:
            del grp['espgrid']
        grp = grp.create_group('espgrid')

        # Store results
        for key, value in results.iteritems():
            grp[key] = value
        ui_grid.to_hdf5(grp.create_group('ui_grid'))


if __name__ == '__main__':
    main()
