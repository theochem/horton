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


# TODO: mere with generic atom script

import sys, argparse, os, stat
from glob import glob
import numpy as np, h5py as h5

from horton import log, lebedev_laikov_npoints, ExpRTransform, RTransform, \
    SimpsonIntegrator1D, AtomicGrid, LebedevLaikovSphereGrid, ProAtomDB, \
    angstrom
from horton.scripts.atomdb import *


def parse_args_input(args):
    parser = argparse.ArgumentParser(prog='horton-atomdb-adf.py input',
        description='Create ADF input files for a database of pro-atoms.',
        epilog=epilog_input)
    add_input_arguments(parser)

    # add specific arguments for the grid
    parser.add_argument(
        "-l", "--lebedev", default=350, type=int,
        help="The number of grid points for the spherical averaging. "
        "[default=%%(default)s]. Select from: %s." % (" ".join(str(i) for i in lebedev_laikov_npoints))
    )
    parser.add_argument(
        "--rmin", default=2e-4, type=float,
        help="The smalest radius for the radial grid (in angstroms). [default=%(default)s]"
    )
    parser.add_argument(
        "--rmax", default=20.0, type=float,
        help="The largest radius for the radial grid (in angstroms). [default=%(default)s]"
    )
    parser.add_argument(
        '-s', "--grid-size", default=100, type=int,
        help="The number of points in the radial grid. [default=%(default)s]"
    )

    return parser.parse_args(args)


grid_prefix = '''\
INPUTFILE TAPE21
UNITS
 Length bohr
END
Grid inline
'''
grid_suffix = '''\
END
density scf
'''


run_adf_script = '''
#!/bin/bash

# make sure adf and the kf2hdf.py script are available before running this
# script.

MISSING=0
if ! which adf &>/dev/null; then echo "adf binary not found."; MISSING=1; fi
if ! which densf &>/dev/null; then echo "densf binary not found."; MISSING=1; fi
if ! which kf2hdf.py &>/dev/null; then echo "kf2hdf.py not found."; MISSING=1; fi
if [ $MISSING -eq 1 ]; then echo "The required programs are not present on your system. Giving up."; exit -1; fi

function do_atom {
    echo "Computing in ${1}"
    cd ${1}
    if [ -e atom.out ]; then
        echo "Output file present in ${1}, not recomputing."
    else
        adf < atom.in > atom.out
        densf < grid.in > grid.out
        kf2hdf.py TAPE41 grid.h5
        rm logfile t21.* TAPE21 TAPE41
    fi
    cd -
}
'''

def main_input(args):
    # Load the template file
    with open(args.template) as f:
        template = Template(f.read())

    # create a grid based on the arguments
    rtf = ExpRTransform(args.rmin*angstrom, args.rmax*angstrom, args.grid_size)
    atspec = (rtf, SimpsonIntegrator1D(), args.lebedev)
    grid = AtomicGrid(0, np.zeros(3), atspec, random_rotate=False)

    # also write the radial grid specification
    with open('rtf.txt', 'w') as f:
        print >> f, rtf.to_string()

    # Loop over all atomic states and make input files
    for number, charge, mult in iter_states(args.elements, args.max_kation, args.max_anion, args.hund):
        dn_mult, do_write = write_input(number, charge, mult, template, args.overwrite)
        if do_write:
            # also write an input file for the grid

            # Write the grid input
            with open('%s/grid.in' % dn_mult, 'w') as f:
                f.write(grid_prefix)
                for point in grid.points:
                    f.write('    %23.16e  %23.16e  %23.16e\n' % tuple(point))
                f.write(grid_suffix)

    # finally write a script
    with open('run_adf.sh', 'w') as f:
        print >> f, run_adf_script
        for dn in sorted(glob("[01]??_??_[01]??_q[+-]??/mult??")):
            print >> f, 'do_atom', dn
    # make the script executable
    os.chmod('run_adf.sh', stat.S_IXUSR | os.stat('run_adf.sh').st_mode)
    if log.do_medium:
        log('Written run_adf.sh')


def parse_args_convert(args):
    parser = argparse.ArgumentParser(prog='horton-atomdb-adf.py convert',
        description='Convert the output of the atomic computations to a horton '
                    'proatom hdf5 file.')
    add_convert_arguments(parser)
    return parser.parse_args(args)


def load_atom_energy(dn_mult):
    with open('%s/atom.out' % dn_mult) as f:
        for line in f:
            if line.startswith('Total Bonding Energy:'):
                return float(line[30:56])


def load_record(dn_mult, rtf):
    with h5.File('%s/grid.h5' % dn_mult) as f:
        npoint = f['Grid/total nr of points'][()]
        assert npoint % rtf.npoint == 0
        nll = npoint / rtf.npoint
        assert nll in lebedev_laikov_npoints
        llgrid = LebedevLaikovSphereGrid(np.zeros(3), 1.0, nll, random_rotate=False)
        dens = f['SCF/Density'][:].reshape((rtf.npoint,nll))
        record = np.dot(dens, llgrid.weights)/llgrid.weights.sum()
        return record


def main_convert(args):
    # load the rtransform
    with open('rtf.txt') as f:
        rtf = RTransform.from_string(f.read().strip())

    # Loop over all sensible directories
    energy_table = EnergyTable()
    records = {}
    for dn_state in sorted(glob("[01]??_??_[01]??_q[+-]??")):
        number = int(dn_state[:3])
        pop = int(dn_state[7:10])

        cases = []
        for dn_mult in sorted(glob('%s/mult??' % dn_state)):
            energy = load_atom_energy(dn_mult)
            if energy is None:
                if log.do_medium:
                    log('No results found:  ', dn_mult)
                continue
            cases.append((energy, dn_mult))
        cases.sort()

        if len(cases) == 0:
            if log.do_medium:
                log('Nothing found in:  ', dn_state)
            continue

        # Get the lowest in energy and write to chk file
        cases.sort()
        energy, dn_mult = cases[0]

        # check for spurious atoms
        if energy_table.is_spurious(number, pop, energy):
            if args.include_spurious:
                if log.do_warning:
                    log.warn('Spurious atom included.')
            else:
                if log.do_medium:
                    log('Skipping spurious: ', dn_state)
                continue


        records[(number, pop)] = load_record(cases[0][1], rtf)
        energy_table.add(number, pop, cases[0][0])
        if log.do_medium:
            log('Succesfull:        ', dn_state)

    if log.do_medium:
        energy_table.log()

    proatomdb = ProAtomDB(rtf, records)
    proatomdb.to_file('atoms.h5')
    if log.do_medium:
        log('Written atoms.h5')

def main():
    args = sys.argv[1:]
    if len(args) == 0:
        print >> sys.stderr, 'Expecting at least one argument: "input" or convert"'
        sys.exit(-1)
    command = args.pop(0)
    if command == 'input':
        parsed = parse_args_input(args)
        main_input(parsed)
    elif command == 'convert':
        parsed = parse_args_convert(args)
        main_convert(parsed)
    else:
        print >> sys.stderr, 'The first argument must be "input" or "convert"'
        sys.exit(-1)


if __name__ == '__main__':
    main()
