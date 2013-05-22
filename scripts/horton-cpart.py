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

from horton import System, cpart_schemes, Cell, ProAtomDB, log, symmetry_analysis
from horton.scripts.common import reduce_data, store_args, parse_pbc, \
    iter_elements, safe_open_h5, write_part_output


def parse_args():
    parser = argparse.ArgumentParser(prog='horton-cpart.py',
        description='Partition the density from a cube file.')

    parser.add_argument('cube',
        help='The cube file.')
    parser.add_argument('scheme', choices=sorted(cpart_schemes.keys()),
        help='The scheme to be used for the partitioning')
    parser.add_argument('atoms',
        help='An HDF5 file with atomic reference densities.')
    parser.add_argument('--stride', default=1, type=int,
        help='Reduce the grid by subsamping with the given stride in all three '
             'directions. Zero and negative values are ignored. '
             '[default=%(default)s]')
    parser.add_argument('--chop', default=0, type=int,
        help='The number of layers to chop of the end of the grid in each '
             'direction. For most codes this should be zero. For Crystal, this '
             'should be 1. [default=%(default)s]')
    parser.add_argument('--overwrite', default=False, action='store_true',
        help='Overwrite existing output in the HDF5 file')
    parser.add_argument('--debug', default=False, action='store_true',
        help='Add additional internal results to a debug subgroup.')
    parser.add_argument('--suffix', default=None, type=str,
        help='Add an additional suffix to the HDF5 output group.')
    parser.add_argument('--pbc', default='111', type=str,
        help='Specify the periodicity. The three digits refer to a, b and c '
             'cell vectors. 1=periodic, 0=aperiodic.')
    parser.add_argument('--symmetry', default=None, type=str,
        help='Perform a symmetry analysis on the AIM results. This option '
             'requires one argument: a CIF file with the generators of the '
             'symmetry of this system and a primitive unit cell.')

    parser.add_argument('--compact', default=None, type=float,
        help='Reduce the cutoff radius of the proatoms such that the tail with '
             'the given number of electrons is neglected. The purpose of this '
             'option is to improve the computational efficiency with a minimal '
             'effect on the results. A typical value is 0.01.')

    parser.add_argument('--wcor', default='1-118', type=str,
        help='The elements for which weight corrections are used. This can be '
             'a comma-separated list of element symbols and/or numbers that '
             'includes ranges. For example, "B,7-9" corresponds to boron, '
             'nitrogen, oxygen and fluorine. The argument 0 will disable '
             'weight corrections entirely.')
    parser.add_argument('--wcor-rcut-max', default=2.0, type=float,
        help='Maximum cutoff radious for weight corrections in Bohr')
    parser.add_argument('--wcor-rcond', default=0.1, type=float,
        help='The regularization strength used for the weight corrections')
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
    fn_h5 = args.cube + '.h5'
    grp_name = '%s_r%i' % (args.scheme, args.stride)
    if args.suffix is not None:
        grp_name += '_'+args.suffix

    if os.path.isfile(fn_h5):
        with safe_open_h5(fn_h5, 'r') as f:
            if 'cpart/%s' % grp_name in f and not args.overwrite:
                if log.do_warning:
                    log.warn('Skipping because "%s" is already present in the output.' % grp_name)
                return

    # Load the system
    sys = System.from_file(args.cube)
    ugrid = sys.props['ugrid']
    ugrid.pbc[:] = parse_pbc(args.pbc)
    moldens = sys.props['cube_data']

    # Reduce the grid if required
    if args.stride > 1 or args.chop > 0:
        moldens, ugrid = reduce_data(moldens, ugrid, args.stride, args.chop)

    # Load the proatomdb and make pro-atoms more compact if that is requested
    proatomdb = ProAtomDB.from_file(args.atoms)
    if args.compact is not None:
        proatomdb.compact(args.compact)
    proatomdb.normalize()

    # Select the partitioning scheme
    CPartClass = cpart_schemes[args.scheme]

    # List of element numbers for which weight corrections are needed:
    wcor_numbers = list(iter_elements(args.wcor))

    # Run the partitioning
    kwargs = dict((key, val) for key, val in vars(args).iteritems() if key in CPartClass.options)
    cpart = cpart_schemes[args.scheme](
        sys, ugrid, True, moldens, proatomdb, wcor_numbers,
        args.wcor_rcut_max, args.wcor_rcond, **kwargs)
    names = cpart.do_all()

    # Do a symmetry analysis if requested.
    if args.symmetry is not None:
        sys_sym = System.from_file(args.symmetry)
        sym = sys_sym.props.get('symmetry')
        if sym is None:
            raise ValueError('No symmetry information found in %s.' % args.symmetry)
        sys_results = dict((name, cpart[name]) for name in names)
        sym_results = symmetry_analysis(sys, sym, sys_results)
        cpart.cache.dump('symmetry', sym_results)
        names.append('symmetry')
        sys.props['symmetry'] = sym

    write_part_output(fn_h5, 'cpart', cpart, grp_name, names, args)

if __name__ == '__main__':
    main()
