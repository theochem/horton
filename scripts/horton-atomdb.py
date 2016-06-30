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


import sys
import argparse
import numpy as np
from glob import glob

from horton import __version__
from horton.log import log
from horton.grid.atgrid import AtomicGridSpec
from horton.part.proatomdb import ProAtomDB, ProAtomRecord
from horton.scripts.atomdb import iter_states, plot_atoms, Template, EnergyTable, \
    atom_programs


# All, except underflows, is *not* fine.
np.seterr(divide='raise', over='raise', invalid='raise')


epilog_input = '''\
The following fields must be present in the template file: ${charge}, ${mult}
and ${element} or ${number}. It may also contain include fields like
${file:some.file} and ${line:some.file}. In the first case, the field will be
replaced by the contents of the file "some.file.NNN_PPP_MM" where NNN is the
element number, PPP is the population and MM is the multiplicity. These numbers
are left-padded with zeros to fix the the length. If a number is zero, it acts
as a wildcard. For example, an include file for any oxygen atom would have
suffix 008_000_00. In the second case, the field is replaced by a single line
from the file "some.file". Each line in this file has a prefix "NNN_PPP_MM ".
The same matching rules are used to select one line from this file. An include
filename may only contain characters from the set [_a-z0-9.-].
'''


def parse_args_input(args):
    parser = argparse.ArgumentParser(prog='horton-atomdb.py input',
        description='Create input files for a database of pro-atoms.',
        epilog=epilog_input)
    parser.add_argument('-V', '--version', action='version',
        version="%%(prog)s (HORTON version %s)" % __version__)

    parser.add_argument('program', choices=sorted(atom_programs.keys()),
        help='The name of the program for which input files must be created.')
    parser.add_argument('elements',
        help='The comma-separated list of elements to be computed. One may '
             'also define ranges with the dash sign. No white space is '
             'allowed. For example, H,6-8 will select hydrogen, carbon, '
             'nitrogen and oxygen.')
    parser.add_argument('template',
        help='A template input file for a single atom in the origin. The '
            'template must contain the fields described below.')
    parser.add_argument('--max-cation', type=int, default=3,
        help='The most positive cation to consider. [default=%(default)s]')
    parser.add_argument('--max-anion', type=int, default=2,
        help='The most negative anion to consider. [default=%(default)s]')
    parser.add_argument('--no-hund', dest='hund', default=True, action='store_false',
        help='When this flag is used, the script does not rely on Hund\'s '
             'rule. Several reasonable multiplicities are tested for each '
             'atom-charge combination and the lowest in energy is selected.')
    parser.add_argument('--overwrite', default=False, action='store_true',
        help='Overwrite existing input files.')

    parsed = parser.parse_args(args)
    parsed.program = atom_programs[parsed.program]
    return parsed


def main_input(args):
    # Load the template file
    with open(args.template) as f:
        template = Template(f.read())

    # Loop over all atomic states and make input files
    for number, charge, mult in iter_states(args.elements, args.max_cation, args.max_anion, args.hund):
        args.program.write_input(number, charge, mult, template, args.overwrite)

    # Write a script that will run all computations (also those created previously)
    args.program.write_run_script()


def parse_args_convert(args):
    parser = argparse.ArgumentParser(prog='horton-atomdb.py convert',
        description='Convert the output of the atomic computations to HORTON '
                    'h5 files.')
    parser.add_argument('-V', '--version', action='version',
        version="%%(prog)s (horton version %s)" % __version__)

    parser.add_argument('--grid', type=str, default='veryfine',
        help='Specify the atomic grid used to construct spherical averages. '
             'Six built-in pruned grids are available: coarse, medium, fine, '
             'veryfine, ultrafine, insane. [default=%(default)s] Not all '
             'elements are supported for each grid type. See documentation for '
             'more details and other possible arguments for this option that '
             'allow a more fine-grained control of the atomic integration '
             'grid. Note that the radial part of this grid is also used for '
             'interpolation in horton-wpart.py')

    return parser.parse_args(args)


def main_convert(args):
    # The atomic grid specification
    agspec = AtomicGridSpec(args.grid)

    # The program is detected based on the run script that is present
    run_scripts = glob("run_*.sh")
    if len(run_scripts) != 1:
        raise RuntimeError('Found %i run_*.sh scripts while exactly one is needed to know which program was used to run the atomic computations.' % len(run_scripts))
    program = atom_programs[run_scripts[0][4:-3]]

    # Loop over all sensible directories
    energy_table = EnergyTable()
    records = []
    for dn_state in sorted(glob("[01]??_??_[01]??_q[+-]??")):
        number = int(dn_state[:3])
        pop = int(dn_state[7:10])

        cases = []
        for dn_mult in sorted(glob('%s/mult??' % dn_state)):
            if log.do_medium:
                log('Loading from', dn_mult)
            data, energy = program.load_atom(dn_mult)
            if energy is None:
                if log.do_medium:
                    log('No (sensible) results found:  ', dn_mult)
                continue
            cases.append((energy, data))

        if len(cases) == 0:
            if log.do_medium:
                log('Nothing found in:  ', dn_state)
            continue

        # Get the lowest in energy and write to chk file
        cases.sort()
        energy, data = cases[0]

        # Add case to energy table
        energy_table.add(number, pop, energy)

        # Write atom to HORTON file if possible
        if data is not None:
            data.to_file('%s/horton.h5' % dn_state)

        # Construct a record for the proatomdb
        records.append(ProAtomRecord.from_iodata(data, agspec))

        # Release memory
        data = None
        del cases

        # Let user know we are alive.
        if log.do_medium:
            log('Succesfull:        ', dn_state)

    # Report energies
    if log.do_medium:
        energy_table.log()

    # Write out atoms file
    proatomdb = ProAtomDB(records)
    proatomdb.to_file('atoms.h5')
    if log.do_medium:
        log('Written atoms.h5')

    # Make nice figures
    plot_atoms(proatomdb)


def main():
    # TODO: work with sub-commands support of argparse.

    args = sys.argv[1:]
    if len(args) == 0:
        print >> sys.stderr, 'Expecting at least one argument: "input" or "convert"'
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
