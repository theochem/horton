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


import sys, argparse, os, numpy as np

from horton import System, wpart_schemes, Cell, ProAtomDB, log, BeckeMolGrid, \
    lebedev_laikov_npoints, AtomicGridSpec, __version__
from horton.scripts.common import get_output_filename, store_args, \
    safe_open_h5, write_part_output


# All, except underflows, is *not* fine.
np.seterr(divide='raise', over='raise', invalid='raise')


def parse_args():
    parser = argparse.ArgumentParser(prog='horton-wpart.py',
        description='Partition the density from a wavefunction file.')
    parser.add_argument('-V', '--version', action='version',
        version="%%(prog)s (horton version %s)" % __version__)

    parser.add_argument('wfn',
        help='The wfn file. Supported formats: fchk, mkl, molden.input')
    parser.add_argument('scheme', choices=sorted(wpart_schemes.keys()),
        help='The scheme to be used for the partitioning')
    parser.add_argument('atoms', default=None, nargs='?',
        help='An HDF5 file with atomic reference densities.')

    parser.add_argument('-o', '--output', default=None,
        help='The HDF5 output filename. The default is to remove the extension '
             'from the input file (if any) and to append \'_wpart.h5\'.')
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
    parser.add_argument('--greedy', default=False, action='store_true',
        help='Keep more precomputed results in memory. This speeds up the '
             'partitioning but consumes more memory. It is only applicable to '
             'the Hirshfeld-I (hi) and Hirhfeld-E (he) schemes.')

    return parser.parse_args()


def main():
    args = parse_args()

    # check if the folder is already present in the output file
    fn_h5 = get_output_filename(args.wfn, 'wpart', args.output)
    grp_name = args.scheme
    if args.suffix is not None:
        grp_name += '_' + args.suffix

    if os.path.isfile(fn_h5):
        with safe_open_h5(fn_h5, 'r') as f:
            if 'wpart/%s' % grp_name in f and not args.overwrite:
                if log.do_warning:
                    log.warn('Skipping because "%s" is already present in the output.' % grp_name)
                return

    # Load the system
    sys = System.from_file(args.wfn)

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
    molgrid = BeckeMolGrid(sys, agspec, mode='only')
    sys.update_grid(molgrid) # for the grid to be written to the output
    wpart = wpart_schemes[args.scheme](sys, molgrid, **kwargs)
    names = wpart.do_all()

    write_part_output(fn_h5, 'wpart', wpart, grp_name, names, args)


if __name__ == '__main__':
    main()
