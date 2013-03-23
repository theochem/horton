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
from horton import System, ESPCost, log


def parse_args():
    parser = argparse.ArgumentParser(prog='horton-esp-fit.py',
        description='Estimate charges from an ESP cost function.')

    parser.add_argument('h5',
        help='[HDF5 filename]:[HDF5 group] The HDF5 file can be created with '
             'horton-esp-cost.py. The group must contain matrices A, B and C '
             'that define the quadratic cost function. The HDF5 file must also '
             'contain a system group in which the atoms are defined with an '
             'array numbers (for the atom numbers) and an array coordinates.')
    parser.add_argument('group',
        help='All results will be stored in this subgroup.')
    parser.add_argument('--qtot', '-q', default=0.0, type=float,
        help='The total charge of the system.')
    parser.add_argument('--ridge', default=0.0, type=float,
        help='The thikonov regularization strength used when solving the '
             'charges.')

    # TODO: more constraint and restraint options
    # TODO: When no group is given in h5 argument, run over all groups that
    #       contain A, B and C.

    return parser.parse_args()


def main():
    args = parse_args()

    # Load the cost function from the HDF5 file
    if args.h5.count(':') != 1:
        raise VallueError('The h5 argument must contain one colon.')
    fn_h5, grp_name = args.h5.split(':')
    with h5.File(fn_h5, 'r') as f:
        grp = f[grp_name]
        sys = System(f['system/coordinates'][:], f['system/numbers'][:])
        cost = ESPCost(grp['A'][:], grp['B'][:], grp['C'][()], sys.natom)

    # Store parameters in output
    results = {}
    results['qtot'] = args.qtot
    results['ridge'] = args.ridge

    # Find the optimal charges
    results['x'] = cost.solve(args.qtot, args.ridge)
    results['charges'] = results['x'][:cost.natom]

    # Related properties
    results['cost'] = cost.value(results['x'])
    if results['cost'] < 0:
        results['rmsd'] = 0.0
    else:
        results['rmsd'] = results['cost']**0.5

    # Worst case stuff
    results['cost_worst'] = cost.worst(args.qtot)
    if results['cost_worst'] < 0:
        results['rmsd_worst'] = 0.0
    else:
        results['rmsd_worst'] = results['cost_worst']**0.5

    # Write some things on screen
    if log.do_medium:
        log('RMSD charges:                  %10.5e' % np.sqrt((results['charges']**2).mean()))
        log('RMSD ESP:                      %10.5e' % results['rmsd'])
        log('Worst RMSD ESP:                %10.5e' % results['rmsd_worst'])
        log.hline()

    # Store the results in an HDF5 file
    with h5.File(fn_h5) as f:
        # Get the group for the output
        grp = f[grp_name]
        if args.group in grp:
            del grp[args.group]
        grp = grp.create_group(args.group)

        # Store results
        for key, value in results.iteritems():
            grp[key] = value


if __name__ == '__main__':
    main()
