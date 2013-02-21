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


from horton import log, lebedev_laikov_npoints, ExpRTransform, \
    RTransform, SimpsonIntegrator1D, AtomicGrid, ProAtomDB, angstrom, periodic
from horton.scripts.atomdb import *
import sys, argparse, numpy as np
from glob import glob


epilog_input = '''\
The following fields must be present in the template file: ${charge}, ${mult}
and ${element} or ${number}. It may also contain include fields like
${file:some.file} and ${line:some.file}. In the first case, the field will be by
the contents of the file "some.file.NNN_PPP_MM" where NNN is the element number,
PPP is the population and MM is the multiplicity. These numbers are left-padded
with zeros to fix the the length. If a number is zero, it acts as a wildcard.
For example, an include file for any oxygen atom would have suffix 008_000_00.
In the second case, the field is replaced by a single line from the file
"some.file". Each line in this file has a prefix "NNN_PPP_MM ". The same
matching rules are used to select one line from this file. An include filename
may only contain characters from the set [_a-z0-9.-].
'''


def parse_args_input(args):
    parser = argparse.ArgumentParser(prog='horton-atomdb.py input',
        description='Create input files for a database of pro-atoms.',
        epilog=epilog_input)

    parser.add_argument('program', choices=atom_programs.keys(),
        help='The name of the program for which input files must be created.')
    parser.add_argument('elements',
        help='The comma-separated list of elements to be computed. One may '
             'also define ranges with the dash sign. No white space is '
             'allowed. For example, 1,6-8 will select hydrogen, carbon, '
             'nitrogen and oxygen.')
    parser.add_argument('template',
        help='A template input file for a single atom in the origin. The '
            'template must contain the fields described below.')
    parser.add_argument('--max-kation', type=int, default=3,
        help='The most positive kation to consider. [default=%(default)s]')
    parser.add_argument('--max-anion', type=int, default=2,
        help='The most negative anion to consider. [default=%(default)s]')
    parser.add_argument('--no-hund', dest='hund', default=True, action='store_false',
        help='When this flag is used, the script does not rely on Hund\'s '
             'rule. Several reasonable multiplicities are tested for each '
             'atom-charge combination and the lowest in energy is selected.')
    parser.add_argument('--overwrite', default=False, action='store_true',
        help='Overwrite existing input files.')
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

    parsed = parser.parse_args(args)
    parsed.program = atom_programs[parsed.program]
    return parsed


def main_input(args):
    # Load the template file
    with open(args.template) as f:
        template = Template(f.read())

    # create a grid based on the arguments and write certain arguments to file
    # for convert call
    rtf = ExpRTransform(args.rmin*angstrom, args.rmax*angstrom, args.grid_size)
    atspec = (rtf, SimpsonIntegrator1D(), args.lebedev)
    atgrid = AtomicGrid(0, np.zeros(3), atspec, random_rotate=False)
    with open('settings.txt', 'w') as f:
        print >> f, rtf.to_string()
        print >> f, args.lebedev
        print >> f, args.program.name

    # Loop over all atomic states and make input files
    for number, charge, mult in iter_states(args.elements, args.max_kation, args.max_anion, args.hund):
        args.program.write_input(number, charge, mult, template, args.overwrite, atgrid)

    # Write a script that will run all computations (also those created previously)
    args.program.write_run_script()


def parse_args_convert(args):
    parser = argparse.ArgumentParser(prog='horton-atomdb.py convert',
        description='Convert the output of the atomic computations to horton '
                    'h5 files.')

    return parser.parse_args(args)


def get_color(index):
    colors = ["#FF0000", "#FFAA00", "#00AA00", "#00AAFF", "#0000FF", "#FF00FF", "#777777"]
    return colors[index%len(colors)]


def plot_atoms(proatomdb):
    try:
        import matplotlib.pyplot as pt
    except ImportError:
        if log.do_warning:
            log.warn('Skipping plots because matplotlib was not found.')
        return

    lss = {True: '-', False: ':'}
    for number in proatomdb.get_numbers():
        r = proatomdb.get_rtransform(number).get_radii()
        symbol = periodic[number].symbol
        charges = proatomdb.get_charges(number)

        # The density (rho)
        pt.clf()
        for i, charge in enumerate(charges):
            record = proatomdb.get_record(number, charge)
            y = record.rho
            ls = lss[record.safe]
            color = get_color(i)
            label = 'q=%+i' % charge
            pt.semilogy(r/angstrom, y, lw=2, ls=ls, label=label, color=color)
        pt.xlim(0, 3)
        pt.ylim(ymin=1e-5)
        pt.xlabel('Distance from the nucleus [A]')
        pt.ylabel('Spherically averaged density [Bohr**-3]')
        pt.title('Proatoms for element %s (%i)' % (symbol, number))
        pt.legend(loc=0)
        fn_png  = 'dens_%03i_%s.png' % (number, symbol)
        pt.savefig(fn_png)
        if log.do_medium:
            log('Written', fn_png)

        # 4*pi*r**2*rho
        pt.clf()
        for i, charge in enumerate(charges):
            record = proatomdb.get_record(number, charge)
            y = record.rho
            ls = lss[record.safe]
            color = get_color(i)
            label = 'q=%+i' % charge
            pt.plot(r/angstrom, 4*np.pi*r**2*y, lw=2, ls=ls, label=label, color=color)
        pt.xlim(0, 3)
        pt.ylim(ymin=0.0)
        pt.xlabel('Distance from the nucleus [A]')
        pt.ylabel('4*pi*r**2*density [Bohr**-1]')
        pt.title('Proatoms for element %s (%i)' % (symbol, number))
        pt.legend(loc=0)
        fn_png  = 'rdens_%03i_%s.png' % (number, symbol)
        pt.savefig(fn_png)
        if log.do_medium:
            log('Written', fn_png)

        # The Fukui functions
        pt.clf()
        for i, charge in enumerate(charges[1:]):
            record0 = proatomdb.get_record(number, charge)
            record1 = proatomdb.get_record(number, charges[i])
            f = record0.rho - record1.rho
            ls = lss[record0.safe and record1.safe]
            color = get_color(i)
            label = '%+i_%+i' % (charge, charges[i])
            pt.semilogy(r/angstrom, f, lw=2, ls=ls, label=label, color=color, alpha=1.0)
            pt.semilogy(r/angstrom, -f, lw=2, ls=ls, color=color, alpha=0.2)
        pt.xlim(0, 3)
        pt.ylim(ymin=1e-5)
        pt.xlabel('Distance from the nucleus [A]')
        pt.ylabel('Fukui function [Bohr**-3]')
        pt.title('Proatoms for element %s (%i)' % (symbol, number))
        pt.legend(loc=0)
        fn_png  = 'fukui_%03i_%s.png' % (number, symbol)
        pt.savefig(fn_png)
        if log.do_medium:
            log('Written', fn_png)

        # 4*pi*r**2*Fukui
        pt.clf()
        for i, charge in enumerate(charges[1:]):
            record0 = proatomdb.get_record(number, charge)
            record1 = proatomdb.get_record(number, charges[i])
            f = record0.rho - record1.rho
            ls = lss[record0.safe and record1.safe]
            color = get_color(i)
            label = '%+i-%+i' % (charge, charges[i])
            pt.plot(r/angstrom, 4*np.pi*r**2*f, lw=2, ls=ls, label=label, color=color)
        pt.xlim(0, 3)
        pt.xlabel('Distance from the nucleus [A]')
        pt.ylabel('4*pi*r**2*Fukui [Bohr**-1]')
        pt.title('Proatoms for element %s (%i)' % (symbol, number))
        pt.legend(loc=0)
        fn_png  = 'rfukui_%03i_%s.png' % (number, symbol)
        pt.savefig(fn_png)
        if log.do_medium:
            log('Written', fn_png)


def main_convert(args):
    # load relevant arguments from input program
    with open('settings.txt') as f:
        rtf = RTransform.from_string(f.next().strip())
        nll = int(f.next())
        program = atom_programs[f.next().strip()]
    int1d = SimpsonIntegrator1D()
    atgrid = AtomicGrid(0, np.zeros(3, float), (rtf, int1d, nll), random_rotate=False, keep_subgrids=1)

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
            system, energy = program.load_atom(dn_mult)
            if energy is None:
                if log.do_medium:
                    log('No results found:  ', dn_mult)
                continue
            cases.append((energy, system, dn_mult))

        if len(cases) == 0:
            if log.do_medium:
                log('Nothing found in:  ', dn_state)
            continue

        # Get the lowest in energy and write to chk file
        cases.sort()
        energy, system, dn_mult = cases[0]

        # Add case to energy table
        energy_table.add(number, pop, energy)

        # Write system to Horton file if possible
        if system is not None:
            system.assign_chk('%s/horton.h5' % dn_state)

        # Construct a record for the proatomdb
        records.append(program.create_record(system, dn_mult, atgrid))

        # Release memory and close h5 files
        system = None
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
