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

import h5py as h5
from horton import System, cpart_schemes, Cell, ProAtomDB, log, Scratch
from horton.scripts.common import reduce_data


def parse_args():
    parser = argparse.ArgumentParser(prog='horton-copart.py',
        description='Partition the density in a cube file.')

    parser.add_argument('cube',
        help='The cube file.')
    parser.add_argument('scheme', choices=cpart_schemes.keys(),
        help='The scheme to be used for the partitioning')
    parser.add_argument('atoms',
        help='An HDF5 file with atomic reference densities.')
    parser.add_argument('--reduce', '-r', default=1, type=int,
        help='Reduce the grid by subsamping with the given stride in all three '
             'directions. Zero and negative values are ignored.')
    parser.add_argument('--overwrite', default=False, action='store_true',
        help='Overwrite existing output in the HDF5 file')
    parser.add_argument('--tmp', type=str, default='.',
        help='A directory where the temporary scratch file can be written.')
    parser.add_argument('--debug', default=False, action='store_true',
        help='Add additional internal results to a debug subgroup.')
    parser.add_argument('--suffix', default=None, type=str,
        help='Add an additional suffix to the HDF5 output group.')

    parser.add_argument('--smooth', '-s', default=False, action='store_true',
        help='Use this option when no special measures are needed to integrate '
             'the cusps accurately.')
    parser.add_argument('--max_iter', '-i', default=100, type=int,
        help='The maximum allowed number of iterations.')
    parser.add_argument('--threshold', '-t', default=1e-4, type=float,
        help='The iterative schemes is converged when the maximum change of '
             'the charges between two iterations drops below this threshold.')

    # TODO: add argument to chop of last slice(s)

    return parser.parse_args()


def main():
    args = parse_args()

    # check if the folder is already present in the output file
    fn_h5 = args.cube + '.h5'
    grp_name = '%s_r%i' % (args.scheme, args.reduce)
    if args.smooth:
        grp_name += '_s'
    if args.suffix is not None:
        grp_name += '_'+args.suffix

    if os.path.isfile(fn_h5):
        with h5.File(fn_h5, 'r') as f:
            if 'cpart/%s' % grp_name in f and not args.overwrite:
                if log.do_warning:
                    log.warn('Skipping because "%s" is already present in the output.' % grp_name)
                return

    # Load the system
    sys = System.from_file(args.cube)
    ui_grid = sys.props['ui_grid']
    moldens = sys.props['cube_data']

    # Reduce the grid if required
    if args.reduce > 1:
        moldens, ui_grid = reduce_data(moldens, ui_grid, args.reduce)

    # Load the proatomdb
    proatomdb = ProAtomDB.from_file(args.atoms)

    # Run the partitioning
    with Scratch('%s/_scratch-PID-%i.h5' % (args.tmp, os.getpid())) as scratch:
        CPartClass = cpart_schemes[args.scheme]
        kwargs = dict((key, val) for key, val in vars(args).iteritems() if key in CPartClass.options)
        cpart = cpart_schemes[args.scheme](sys, ui_grid, moldens, proatomdb, scratch, **kwargs)
        names = cpart.do_all()

    # Store the results in an HDF5 file
    with h5.File(fn_h5) as f:
        # Store essential system info
        # TODO: this should be implemented with an improved implementation of System.to_file
        sys_grp = f.require_group('system')
        coordinates = sys_grp.require_dataset('coordinates', sys.coordinates.shape, float, exact=True)
        coordinates[:] = sys.coordinates
        numbers = sys_grp.require_dataset('numbers', sys.numbers.shape, long, exact=True)
        numbers[:] = sys.numbers
        if 'cell' in sys_grp:
            del sys_grp['cell']
        cell_grp = sys_grp.create_group('cell')
        sys.cell.to_hdf5(cell_grp)

        # Store results
        grp_cpart = f.require_group('cpart')
        if grp_name in grp_cpart:
            del grp_cpart[grp_name]
        grp = grp_cpart.create_group(grp_name)
        grp['grid_rvecs'] = ui_grid.grid_cell.rvecs
        for name in names:
            grp[name] = cpart[name]

        if args.debug:
            # Store additional data for debugging
            if 'debug' in grp:
                del grp['debug']
            grp_debug = grp.create_group('debug')
            for debug_key in cpart._cache._store:
                debug_name = '_'.join(str(x) for x in debug_key)
                if debug_name not in names:
                    grp_debug[debug_name] = cpart._cache.load(*debug_key)


if __name__ == '__main__':
    main()
