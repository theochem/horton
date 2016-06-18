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


import argparse, os, numpy as np

from horton import IOData, Cell, ProAtomDB, log, BeckeMolGrid, \
    lebedev_laikov_npoints, AtomicGridSpec, __version__
from horton.scripts.common import store_args, write_part_output, parse_h5, \
    check_output
from horton.scripts.wpart import wpart_slow_analysis, wpart_schemes


# All, except underflows, is *not* fine.
np.seterr(divide='raise', over='raise', invalid='raise')


def parse_args():
    parser = argparse.ArgumentParser(prog='horton-wpart.py',
        description='Partition the density from a wavefunction file.')
    parser.add_argument('-V', '--version', action='version',
        version="%%(prog)s (HORTON version %s)" % __version__)

    parser.add_argument('wfn',
        help='The wfn file. Supported formats: fchk, mkl, molden.input, wfn')
    parser.add_argument('output',
        help='The output destination in the form file.h5:group. The colon and '
             'the group name are optional. When omitted, the root group of the '
             'HDF5 file is used. '
             'To mimick the behavior of HORTON 1.2.0 scripts, use '
             '"${prefix}_wpart.h5:wpart/${scheme}" where ${prefix} '
             'is the name of the wfn file without extension and ${scheme} is '
             'the partitioning scheme used.')
    parser.add_argument('scheme', choices=sorted(wpart_schemes.keys()),
        help='The scheme to be used for the partitioning')
    parser.add_argument('atoms', default=None, nargs='?',
        help='An HDF5 file with atomic reference densities.')

    parser.add_argument('--overwrite', default=False, action='store_true',
        help='Overwrite existing output in the HDF5 file')
    parser.add_argument('--debug', default=False, action='store_true',
        help='Add additional internal results to a debug subgroup.')
    parser.add_argument('--suffix', default=None, type=str,
        help='Add an additional suffix to the HDF5 output group.')

    parser.add_argument('--grid', type=str, default='medium',
        help='Specify the atomic integration grids. Six built-in pruned '
             'grids are available: coarse, medium, fine, veryfine, ultrafine, '
             'insane. [default=%(default)s] Not all elements are supported '
             'for each grid type. See documentation for more details and other '
             'possible arguments for this option that allow a more '
             'fine-grained control of the atomic integration grid.')
    parser.add_argument('-e', '--epsilon', default=1e-8, type=float,
        help='Allow errors on the computed electron density of this magnitude '
             'for the sake of efficiency.')
    parser.add_argument('--maxiter', '-i', default=500, type=int,
        help='The maximum allowed number of iterations. [default=%(default)s]')
    parser.add_argument('--threshold', '-t', default=1e-6, type=float,
        help='The iterative scheme is converged when the maximum change of '
             'the charges between two iterations drops below this threshold. '
             '[default=%(default)s]')
    parser.add_argument('--lmax', default=3, type=int,
        help='The maximum angular momentum to consider in multipole expansions')
    parser.add_argument('--slow', default=False, action='store_true',
        help='Also compute the more expensive AIM properties that require the '
             'AIM overlap matrices.')

    return parser.parse_args()


def main():
    args = parse_args()

    fn_h5, grp_name = parse_h5(args.output, 'output')
    # check if the group is already present (and not empty) in the output file
    if check_output(fn_h5, grp_name, args.overwrite):
        return

    # Load the system
    mol = IOData.from_file(args.wfn)

    # Define a list of optional arguments for the WPartClass:
    WPartClass = wpart_schemes[args.scheme]
    kwargs = dict((key, val) for key, val in vars(args).iteritems() if key in WPartClass.options)

    # Load the proatomdb
    if args.atoms is not None:
        proatomdb = ProAtomDB.from_file(args.atoms)
        proatomdb.normalize()
        kwargs['proatomdb'] = proatomdb
    else:
        proatomdb = None

    # Run the partitioning
    agspec = AtomicGridSpec(args.grid)
    grid = BeckeMolGrid(mol.coordinates, mol.numbers, mol.pseudo_numbers, agspec, mode='only')
    dm_full = mol.get_dm_full()
    moldens = mol.obasis.compute_grid_density_dm(dm_full, grid.points, epsilon=args.epsilon)
    dm_spin = mol.get_dm_spin()
    if dm_spin is not None:
        kwargs['spindens'] = mol.obasis.compute_grid_density_dm(dm_spin, grid.points, epsilon=args.epsilon)
    wpart = wpart_schemes[args.scheme](mol.coordinates, mol.numbers, mol.pseudo_numbers,grid, moldens, **kwargs)
    keys = wpart.do_all()

    if args.slow:
        # ugly hack for the slow analysis involving the AIM overlap operators.
        wpart_slow_analysis(wpart, mol)
        keys = list(wpart.cache.iterkeys(tags='o'))

    write_part_output(fn_h5, grp_name, wpart, keys, args)


if __name__ == '__main__':
    main()
