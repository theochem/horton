#!/usr/bin/env python
# -*- coding: utf-8 -*-
# HORTON: Helpful Open-source Research TOol for N-fermion systems.
# Copyright (C) 2011-2017 The HORTON Development Team
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


import argparse
import numpy as np

from horton import log, __version__
from horton.scripts.common import parse_h5, write_script_output, \
    check_output
from horton.scripts.espfit import load_cost


# All, except underflows, is *not* fine.
np.seterr(divide='raise', over='raise', invalid='raise')


def parse_args():
    parser = argparse.ArgumentParser(prog='horton-esp-fit.py',
        description='Estimate charges from an ESP cost function.')
    parser.add_argument('-V', '--version', action='version',
        version="%%(prog)s (HORTON version %s)" % __version__)

    parser.add_argument('cost',
        help='The location of the cost function in the form '
             '"file.h5:group/cost". This argument must be the same as the '
             'output argument of the script horton-esp-cost.py.')
    parser.add_argument('output',
        help='The output destination in the form file.h5:group. The colon and '
             'the group name are optional. When omitted, the root group of the '
             'HDF5 file is used.')
    parser.add_argument('--overwrite', default=False, action='store_true',
        help='Overwrite existing output in the HDF5 file')

    parser.add_argument('--qtot', '-q', default=0.0, type=float,
        help='The total charge of the system. [default=%(default)s]')
    parser.add_argument('--ridge', default=0.0, type=float,
        help='The thikonov regularization strength used when solving the '
             'charges. [default=%(default)s]')

    # TODO: more constraint and restraint options

    return parser.parse_args()


def main():
    args = parse_args()

    fn_h5, grp_name = parse_h5(args.output, 'output')
    # check if the group is already present (and not empty) in the output file
    if check_output(fn_h5, grp_name, args.overwrite):
        return

    # Load the cost function from the HDF5 file
    cost, used_volume = load_cost(args.cost)

    # Find the optimal charges
    results = {}
    results['x'] = cost.solve(args.qtot, args.ridge)
    results['charges'] = results['x'][:cost.natom]

    # Related properties
    results['cost'] = cost.value(results['x'])
    if results['cost'] < 0:
        results['rmsd'] = 0.0
    else:
        results['rmsd'] = (results['cost']/used_volume)**0.5

    # Worst case stuff
    results['cost_worst'] = cost.worst(0.0)
    if results['cost_worst'] < 0:
        results['rmsd_worst'] = 0.0
    else:
        results['rmsd_worst'] = (results['cost_worst']/used_volume)**0.5

    # Write some things on screen
    if log.do_medium:
        log('Important parameters:')
        log.hline()
        log('RMSD charges:                  %10.5e' % np.sqrt((results['charges']**2).mean()))
        log('RMSD ESP:                      %10.5e' % results['rmsd'])
        log('Worst RMSD ESP:                %10.5e' % results['rmsd_worst'])
        log.hline()

    # Store the results in an HDF5 file
    write_script_output(fn_h5, grp_name, results, args)


if __name__ == '__main__':
    main()
