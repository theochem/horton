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


from horton import log, System
from horton.scripts.atomdb import *
import sys, argparse, os
from glob import glob


def parse_args_input(args):
    parser = argparse.ArgumentParser(prog='horton-atomdb.py input',
        description='Create input files for a database of pro-atoms.',
        epilog=epilog_input)
    add_input_arguments(parser)
    return parser.parse_args(args)


def main_input(args):
    # Load the template file
    with open(args.template) as f:
        template = Template(f.read())

    # Loop over all atomic states and make input files
    for number, charge, mult in iter_states(args.elements, args.max_kation, args.max_anion, args.hund):
        write_input(number, charge, mult, template, args.overwrite)


def parse_args_convert(args):
    parser = argparse.ArgumentParser(prog='horton-atomdb.py convert',
        description='Convert the output of the atomic computations to horton '
                    'h5 files.')
    add_convert_arguments(parser)
    return parser.parse_args(args)


def get_fn(dn, pattern):
    fns = sorted(glob('%s/%s' % (dn, pattern)))
    if len(fns) == 0:
        return None
    elif len(fns) == 1:
        return fns[0]
    elif log.do_warning:
        log.warn('Found more than one "%s" file in %s. This is suspicious. Skipping directory.' % (pattern, dn))
        return False


def load_atom_system(dn_mult):
    for pattern in '*.fchk', '*.mkl', '*.molden.input':
        fn = get_fn(dn_mult, pattern)
        if fn is False:
            return None
        elif fn is not None:
            try:
                return System.from_file(fn)
            except:
                if log.do_warning:
                    log.warn('Some error occured while loading "%s"' % fn)
                raise


def main_convert(args):
    # Loop over all sensible directories
    energy_table = EnergyTable()
    for dn_state in sorted(glob("[01]??_??_[01]??_q[+-]??")):
        number = int(dn_state[:3])
        pop = int(dn_state[7:10])

        systems = []
        for dn_mult in sorted(glob('%s/mult??' % dn_state)):
            system = load_atom_system(dn_mult)
            if system is None:
                if log.do_medium:
                    log('No results found:  ', dn_mult)
                continue
            systems.append(system)

        if len(systems) == 0:
            if log.do_medium:
                log('Nothing found in:  ', dn_state)
            continue

        # Get the lowest in energy and write to chk file
        systems.sort(key=(lambda s: s.props['energy']))

        assert pop == systems[0].wfn.nel
        energy = systems[0].props['energy']
        # check for spurious atoms
        if energy_table.is_spurious(number, pop, energy):
            if args.include_spurious:
                if log.do_warning:
                    log.warn('Spurious atom included.')
            else:
                if log.do_medium:
                    log('Skipping spurious: ', dn_state)
                continue

        energy_table.add(number, pop, energy)

        if log.do_medium:
            log('Succesfull:        ', dn_state)
        systems[0].assign_chk('%s/horton.h5' % dn_state)
        del systems

    if log.do_medium:
        energy_table.log()


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
