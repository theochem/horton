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
from horton.scripts.common import parse_h5


def parse_args():
    parser = argparse.ArgumentParser(prog='horton-esp-test.py',
        description='Test how well charges reproduce the ESP.')

    parser.add_argument('h5_esp',
        help='[HDF5 filename]:[HDF5 group] The HDF5 file can be created with '
             'horton-esp-cost.py. The group must contain matrices A, B and C '
             'that define the quadratic cost function. The HDF5 file must also '
             'contain a system group in which the atoms are defined with an '
             'array numbers (for the atom numbers) and an array coordinates.')
    parser.add_argument('h5_charges',
        help='[HDF5 filename]:[HDF5 group] The HDF5 file with the charges. '
             'The results of the test will be written in a subgroup "esp". '
             'The HDF5 file must also contain a system group that matches with '
             'the system group of the h5_esp argument.')
    parser.add_argument('--qtot', '-q', default=None, type=float,
        help='The total charge of the system. When given, the charges from the '
             'HDF5 file are corrected.')

    # TODO: When no group is given in an h5 argument, loop over all relevant groups.

    return parser.parse_args()


def main():
    args = parse_args()

    # Load the cost function from the first HDF5 file
    fn_h5, grp_name = parse_h5(args.h5_esp)
    with h5.File(fn_h5, 'r') as f:
        grp = f[grp_name]
        sys = System(f['system/coordinates'][:], f['system/numbers'][:])
        cost = ESPCost(grp['A'][:], grp['B'][:], grp['C'][()], sys.natom)

    # Load the charges from the second HDF5 file
    fn_h5, grp_name = parse_h5(args.h5_charges)
    with h5.File(fn_h5, 'r') as f:
        grp = f[grp_name]
        charges = grp['charges'][:]
        assert (sys.coordinates == f['system/coordinates'][:]).all()
        assert (sys.numbers == f['system/numbers'][:]).all()

    # Fix total charge if requested
    if args.qtot is not None:
        charges -= (charges.sum() - args.qtot)/len(charges)

    # Store parameters in output
    results = {}
    results['qtot'] = charges.sum()

    # Fitness of the charges
    results['cost'] = cost.value_charges(charges)
    if results['cost'] < 0:
        results['rmsd'] = 0.0
    else:
        results['rmsd'] = results['cost']**0.5

    # Worst case stuff
    results['cost_worst'] = cost.worst(charges.sum())
    if results['cost_worst'] < 0:
        results['rmsd_worst'] = 0.0
    else:
        results['rmsd_worst'] = results['cost_worst']**0.5

    # Write some things on screen
    if log.do_medium:
        log('RMSD charges:                  %10.5e' % np.sqrt((charges**2).mean()))
        log('RMSD ESP:                      %10.5e' % results['rmsd'])
        log('Worst RMSD ESP:                %10.5e' % results['rmsd_worst'])
        log.hline()

    # Store the results in an HDF5 file
    with h5.File(fn_h5) as f:
        # Get the group for the output
        grp = f[grp_name]
        if 'esp' in grp:
            del grp['esp']
        grp = grp.create_group('esp')

        # Store results
        for key, value in results.iteritems():
            grp[key] = value


if __name__ == '__main__':
    main()
