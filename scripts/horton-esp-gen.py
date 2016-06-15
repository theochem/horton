#!/usr/bin/env python
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


import sys, argparse, os, numpy as np

from horton import IOData, UniformGrid, log, angstrom, \
    compute_esp_grid_cube, __version__
from horton.scripts.common import parse_h5, parse_ewald_args, store_args, \
    check_output, write_script_output
from horton.scripts.espfit import load_charges


# All, except underflows, is *not* fine.
np.seterr(divide='raise', over='raise', invalid='raise')


def parse_args():
    parser = argparse.ArgumentParser(prog='horton-esp-gen.py',
        description='Generate electrostatic potential grid data from charges '
                    'for a 3D periodic system.')
    parser.add_argument('-V', '--version', action='version',
        version="%%(prog)s (HORTON version %s)" % __version__)

    parser.add_argument('charges', type=str,
        help='The atomic charges to be used in the form '
             '"file.h5:group/charges". ')
    parser.add_argument('grid', type=str,
        help='Any type of file that contains a uniform grid specification and '
             'atomic coordinates, e.g. a Gaussian cube file.')
    parser.add_argument('output', type=str,
        help='The output destination in the form file.h5:group. The colon and '
             'the group name are optional. When omitted, the root group of the '
             'HDF5 file is used.')
    parser.add_argument('--overwrite', default=False, action='store_true',
        help='Overwrite existing output in the HDF5 file')

    parser.add_argument('--qtot', '-q', default=None, type=float,
        help='The total charge of the system. When given, the charges from the '
             'HDF5 file are corrected after they are read in.')
    parser.add_argument('--rcut', default=10.0, type=float,
        help='The real-space cutoff for the electrostatic interactions in '
             'angstrom. [default=%(default)s]')
    parser.add_argument('--alpha-scale', default=3.0, type=float,
        help='The alpha scale (alpha = alpha_scale/rcut) for the separation '
             'between short-range and long-range electrostatic interactions. '
             '[default=%(default)s]')
    parser.add_argument('--gcut-scale', default=1.1, type=float,
        help='The gcut scale (gcut = gcut_scale*alpha) for the reciprocal '
             'space constribution to the electrostatic interactions. '
             '[default=%(default)s]')

    return parser.parse_args()


def load_ugrid_coordinates(arg_grid):
    mol = IOData.from_file(arg_grid)
    return mol.grid, mol.coordinates


def main():
    args = parse_args()

    fn_h5, grp_name = parse_h5(args.output, 'output')
    # check if the group is already present (and not empty) in the output file
    if check_output(fn_h5, grp_name, args.overwrite):
        return

    # Load the charges from the HDF5 file
    charges = load_charges(args.charges)

    # Load the uniform grid and the coordintes
    ugrid, coordinates = load_ugrid_coordinates(args.grid)
    ugrid.pbc[:] = 1 # enforce 3D periodic

    # Fix total charge if requested
    if args.qtot is not None:
        charges -= (charges.sum() - args.qtot)/len(charges)

    # Store parameters in output
    results = {}
    results['qtot'] = charges.sum()

    # Determine the grid specification
    results['ugrid'] = ugrid

    # Ewald parameters
    rcut, alpha, gcut = parse_ewald_args(args)

    # Some screen info
    if log.do_medium:
        log('Important parameters:')
        log.hline()
        log('Number of grid points:   %12i' % ugrid.size)
        log('Grid shape:                 [%8i, %8i, %8i]' % tuple(ugrid.shape))
        log('Ewald real cutoff:       %12.5e' % rcut)
        log('Ewald alpha:             %12.5e' % alpha)
        log('Ewald reciprocal cutoff: %12.5e' % gcut)
        log.hline()
        # TODO: add summation ranges
        log('Computing ESP (may take a while)')

    # Allocate and compute ESP grid
    esp = np.zeros(ugrid.shape, float)
    compute_esp_grid_cube(ugrid, esp, coordinates, charges, rcut, alpha, gcut)
    results['esp'] = esp

    # Store the results in an HDF5 file
    write_script_output(fn_h5, grp_name, results, args)


if __name__ == '__main__':
    main()
