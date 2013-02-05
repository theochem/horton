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


from horton import *
import sys, argparse, os
from glob import glob
from string import Template


log.set_level(log.silent)


def iter_elements(elements_str):
    for item in elements_str.split(','):
        if '-' in item:
            words = item.split("-")
            if len(words) != 2:
                raise ValueError("Each item should contain at most one dash.")
            begin = int(words[0])
            end = int(words[1])
            for number in xrange(begin,end+1):
                yield number
        else:
            yield int(item)


def parse_args_input(args):
    parser = argparse.ArgumentParser(prog='horton-atomdb.py input')
    parser.add_argument('elements',
        help='The comma-separated list of elements to be computed. One may '
             'also define ranges with the dash sign. No white space is '
             'allowed. For example, 1,6-8 will select hydrogen, carbon, '
             'nitrogen and oxygen.')
    parser.add_argument('template',
        help='A template input file for a single atom in the origin. It must '
             'contain fields ${charge}, ${mult} and ${element} or ${number}.')
    parser.add_argument('--max-kation', type=int, default=3,
        help='The most positive kation to consider. [default=3]')
    parser.add_argument('--max-anion', type=int, default=2,
        help='The most negative anion to consider. [default=2]')
    parser.add_argument('--no-hund', dest='hund', default=True, action='store_false',
        help='When this flag is used, the script does not rely on Hund\'s '
             'rule. Several reasonable multiplicities are tested for each '
             'atom-charge combination and the lowest in energy is selected.')
    parser.add_argument('--overwrite', default=False, action='store_true',
        help='Overwrite existing input files.')
    return parser.parse_args(args)


# Presets for spin multiplicites. The first element is according to Hund's rule.
# Following elements are reasonable.
mult_presets = {
    1: [2],
    2: [1, 3],
    3: [2, 4],
    4: [1, 3],
    5: [2, 4],
    6: [3, 5, 1],
    7: [4, 2],
    8: [3, 1],
    9: [2],
    10: [1],
    11: [2],
    12: [1, 3],
    13: [2, 4],
    14: [3, 5, 1],
    15: [4, 2],
    16: [3, 1],
    17: [2],
    18: [1],
    19: [2],
    20: [1, 3],
    21: [2, 4],
    22: [3, 5, 1],
    23: [4, 6, 2],
    24: [7, 5, 3, 1],
    25: [6, 4, 2],
    26: [5, 3, 1],
    27: [4, 2],
    28: [3, 1],
    29: [2],
    30: [1],
    31: [2, 4, 6],
    32: [3, 1, 5],
    33: [4, 2],
    34: [3, 1],
    35: [2],
    36: [1],
    37: [2],
    38: [1, 3],
    39: [2, 4],
    40: [3, 1, 5],
    41: [6, 4, 2],
    42: [7, 5, 3, 1],
    43: [6, 4, 2],
    44: [5, 3, 1],
    45: [4, 2],
    46: [1, 3],
    47: [2],
    48: [1],
    49: [2, 4],
    50: [3, 1, 5],
    51: [4, 2],
    52: [3, 1],
    53: [2],
    54: [1],
    55: [2],
    56: [1, 3],
}


def iter_mults(nel, hund):
    if hund:
        yield mult_presets[nel][0]
    else:
        for mult in mult_presets[nel]:
            yield mult


def main_input(args):
    # Load the template file
    with open(args.template) as f:
        template = Template(f.read())
    base_inp = os.path.basename(args.template)

    # Loop over all elements and make input files
    for number in iter_elements(args.elements):
        for charge in xrange(-args.max_anion, args.max_kation+1):
            nel = number - charge
            if nel <= 0:
                continue
            dn_state = '%03i_%s_%03i_q%+03i' % (
                number, periodic[number].symbol.lower().rjust(2, '_'), nel, charge)
            if not os.path.isdir(dn_state):
                os.mkdir(dn_state)
            for mult in iter_mults(nel, args.hund):
                dn_mult = '%s/mult%02i' % (dn_state, mult)
                if not os.path.isdir(dn_mult):
                    os.mkdir(dn_mult)
                fn_inp = '%s/%s' % (dn_mult, base_inp)
                exists = os.path.isfile(fn_inp)
                if not exists or args.overwrite:
                    with open(fn_inp, 'w') as f:
                        f.write(template.safe_substitute(
                            charge=str(charge),
                            mult=str(mult),
                            number=str(number),
                            element=periodic[number].symbol,
                        ))
                    if exists:
                        print 'Overwritten:      ', fn_inp
                    else:
                        print 'Written new:      ', fn_inp
                else:
                    print 'Not overwriting:  ', fn_inp


def parse_args_convert(args):
    parser = argparse.ArgumentParser(prog='horton-atomdb.py convert')
    return parser.parse_args(args)


def get_fn(dn, pattern):
    fns = sorted(glob('%s/%s' % (dn, pattern)))
    if len(fns) == 0:
        return None
    elif len(fns) == 1:
        return fns[0]
    else:
        print 'Found more than one "%s" file in %s. This is suspicious. Skipping directory.' % (pattern, dn)
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
                print 'Some error occured while loading "%s"' % fn
                raise


def main_convert(args):
    # Loop over all sensible directories
    for dn_state in sorted(glob("[01]??_??_[01]??_q[+-]??")):
        systems = []
        for dn_mult in sorted(glob('%s/mult??' % dn_state)):
            system = load_atom_system(dn_mult)
            if system is None:
                print 'No results found:  ', dn_mult
                continue
            systems.append(system)

        if len(systems) == 0:
            print 'Nothing found in:  ', dn_state
            continue

        # Get the lowest in energy and write to chk file
        systems.sort(key=(lambda s: s.props['energy']))
        print 'Succesfull:        ', dn_state
        systems[0].assign_chk('%s/horton.h5' % dn_state)
        del systems


def main():
    args = sys.argv[1:]
    if len(args) == 0:
        print >> sys.stderr, 'Expecting at least one argument: "input" or convert"'
    command = args.pop(0)
    if command == 'input':
        parsed = parse_args_input(args)
        main_input(parsed)
    elif command == 'convert':
        parsed = parse_args_convert(args)
        main_convert(parsed)
    else:
        print >> sys.stderr, 'The first argument must be "input" or "convert"'


if __name__ == '__main__':
    main()
