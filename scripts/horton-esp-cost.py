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

from horton import IOData, setup_weights, ESPCost, log, \
    angstrom, __version__
from horton.scripts.common import reduce_data, parse_ewald_args, parse_pbc, \
    write_script_output, parse_h5, check_output
from horton.scripts.espfit import parse_wdens, parse_wnear, parse_wfar, \
    load_rho, max_at_edge
from horton.grid.cext import UniformGrid


# All, except underflows, is *not* fine.
np.seterr(divide='raise', over='raise', invalid='raise')


def parse_args():
    parser = argparse.ArgumentParser(prog='horton-esp-cost.py',
        description='Construct a cost function and fit charges to ESP.')
    parser.add_argument('-V', '--version', action='version',
        version="%%(prog)s (HORTON version %s)" % __version__)

    parser.add_argument('cube',
        help='The cube file.')
    parser.add_argument('output',
        help='The output destination in the form file.h5:group. The colon and '
             'the group name are optional. When omitted, the root group of the '
             'HDF5 file is used. '
             'To mimick the behavior of HORTON 1.2.0 scripts, use '
             '"${prefix}_espfit.h5:espfit/cost_r${stride}" where ${prefix} '
             'is the name of the cube file without extension and ${stride} is '
             'the value of the stride option.')
    parser.add_argument('--overwrite', default=False, action='store_true',
        help='Overwrite existing output in the HDF5 file')

    parser.add_argument('--sign', '-s', default=False, action='store_true',
        help='Change the sign of the ESP data. This is needed for CP2K and VASP.')
    parser.add_argument('--stride', default=1, type=int,
        help='Reduce the grid by subsamping with the given stride in all three '
             'directions. Zero and negative values are ignored. '
             '[default=%(default)s]')
    parser.add_argument('--chop', default=0, type=int,
        help='The number of layers to chop of the end of the grid in each '
             'direction. For most codes this should be zero. For Crystal, this '
             'should be 1.')

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

    parser.add_argument('--wdens', default=None, type=str, nargs='?', const=':-9:0.8',
        help='Define weights based on an electron density. The argument has '
             'the following format: "dens.cube:lnrho0:sigma". All arguments '
             'are optional. The density is loaded from the file in the '
             'first field or a prodensity is constructed if the file is not '
             'given. The cube file must have the same header is the potential '
             'data file. The second and third field are parameters of the '
             'wieght function and are defined in DOI:10.1021/ct600295n. '
             'The default is \':-9:0.8\'.')
    parser.add_argument('--wnear', default=None, type=str, nargs='+',
        help='Define weights that go to zero near the nuclei. Multiple '
             'arguments are allowed and have the following format: '
             '"number:r0:gamma". The first field is the elemenet number and '
             'the last two fields define the switching function and are in '
             'Angstrom. The last field is optional and is 0.5*A by default. '
             'The second field is the middle of the switching function. '
             'The third field is the half width of the switching function. If '
             'number is zero, the argument applies to all elements.')
    parser.add_argument('--wfar', default=None, type=str,
        help='Define weights that go to zero far from the nuclei. The '
             'argument has the following format: "r0:gamma". '
             'The first field is the distance (to the closest nucleus) at '
             'which the weight switches to zero. The second field is the half '
             'width of the switching function. To avoid artifacts, a smooth '
             'version of "distance to the closest nucleus" is used. The second '
             'field is optional and has 1.0 as default value')
    parser.add_argument('--wsave', default=None, type=str,
        help='Save the weights array to the given cube file.')
    parser.add_argument('--pbc', default='111', type=str, choices=['000', '111'],
        help='Specify the periodicity. The three digits refer to a, b and c '
             'cell vectors. 1=periodic, 0=aperiodic.')

    return parser.parse_args()


def main():
    args = parse_args()

    fn_h5, grp_name = parse_h5(args.output, 'output')
    # check if the group is already present (and not empty) in the output file
    if check_output(fn_h5, grp_name, args.overwrite):
        return

    # Load the potential data
    if log.do_medium:
        log('Loading potential array')
    mol_pot = IOData.from_file(args.cube)
    if not isinstance(mol_pot.grid, UniformGrid):
        raise TypeError('The specified file does not contain data on a rectangular grid.')
    mol_pot.grid.pbc[:] = parse_pbc(args.pbc) # correct pbc
    esp = mol_pot.cube_data

    # Reduce the grid if required
    if args.stride > 1:
        esp, mol_pot.grid = reduce_data(esp, mol_pot.grid, args.stride, args.chop)

    # Fix sign
    if args.sign:
        esp *= -1

    # Some screen info
    if log.do_medium:
        log('Important parameters:')
        log.hline()
        log('Number of grid points:   %12i' % np.product(mol_pot.grid.shape))
        log('Grid shape:                 [%8i, %8i, %8i]' % tuple(mol_pot.grid.shape))
        log('PBC:                        [%8i, %8i, %8i]' % tuple(mol_pot.grid.pbc))
        log.hline()

    # Construct the weights for the ESP Cost function.
    wdens = parse_wdens(args.wdens)
    if wdens is not None:
        if log.do_medium:
            log('Loading density array')
        # either the provided density or a built-in prodensity
        rho = load_rho(mol_pot.coordinates, mol_pot.numbers, wdens[0], mol_pot.grid, args.stride, args.chop)
        wdens = (rho,) + wdens[1:]
    if log.do_medium:
        log('Constructing weight function')
    weights = setup_weights(mol_pot.coordinates, mol_pot.numbers, mol_pot.grid,
        dens=wdens,
        near=parse_wnear(args.wnear),
        far=parse_wnear(args.wfar),
    )

    # write the weights to a cube file if requested
    if args.wsave is not None:
        if log.do_medium:
            log('   Saving weights array   ')
        # construct a new data dictionary that contains all info for the cube file
        mol_weights = mol_pot.copy()
        mol_weights.cube_data = weights
        mol_weights.to_file(args.wsave)

    # rescale weights such that the cost function is the mean-square-error
    if weights.max() == 0.0:
        raise ValueError('No points with a non-zero weight were found')
    wmax = weights.min()
    wmin = weights.max()
    used_volume = mol_pot.grid.integrate(weights)

    # Some screen info
    if log.do_medium:
        log('Important parameters:')
        log.hline()
        log('Used number of grid points:   %12i' % (weights>0).sum())
        log('Used volume:                      %12.5f' % used_volume)
        log('Used volume/atom:                 %12.5f' % (used_volume/mol_pot.natom))
        log('Lowest weight:                %12.5e' % wmin)
        log('Highest weight:               %12.5e' % wmax)
        log('Max weight at edge:           %12.5f' % max_at_edge(weights, mol_pot.grid.pbc))

    # Ewald parameters
    rcut, alpha, gcut = parse_ewald_args(args)

    # Some screen info
    if log.do_medium:
        log('Ewald real cutoff:       %12.5e' % rcut)
        log('Ewald alpha:             %12.5e' % alpha)
        log('Ewald reciprocal cutoff: %12.5e' % gcut)
        log.hline()

    # Construct the cost function
    if log.do_medium:
        log('Setting up cost function (may take a while)   ')
    cost = ESPCost.from_grid_data(mol_pot.coordinates, mol_pot.grid, esp, weights, rcut, alpha, gcut)

    # Store cost function info
    results = {}
    results['cost'] = cost
    results['used_volume'] = used_volume

    # Store cost function properties
    results['evals'] = np.linalg.eigvalsh(cost._A)
    abs_evals = abs(results['evals'])
    if abs_evals.min() == 0.0:
        results['cn'] = 0.0
    else:
        results['cn'] = abs_evals.max()/abs_evals.min()

    # Report some on-screen info
    if log.do_medium:
        log('Important parameters:')
        log.hline()
        log('Lowest abs eigen value:       %12.5e' % abs_evals.min())
        log('Highest abs eigen value:      %12.5e' % abs_evals.max())
        log('Condition number:             %12.5e' % results['cn'])
        log.hline()

    # Store the results in an HDF5 file
    write_script_output(fn_h5, grp_name, results, args)


if __name__ == '__main__':
    main()
