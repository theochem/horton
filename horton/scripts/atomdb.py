# -*- coding: utf-8 -*-
# Horton is a Density Functional Theory program.
# Copyright (C) 2011-2013 Toon Verstraelen <Toon.Verstraelen@UGent.be>
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


from string import Template as BaseTemplate
import re, os

from horton import periodic, log


__all__ = [
    'iter_elements', 'iter_mults', 'Template', 'add_input_arguments',
    'epilog_input', 'add_convert_arguments', 'iter_states', 'write_input',
    'EnergyTable',
]


def iter_elements(elements_str):
    '''Interpret a string as a list of elements

       elements_str
            A string with comma-separated element numbers. One may add ranges
            with the format 'N-M' where M>N.
    '''
    for item in elements_str.split(','):
        if '-' in item:
            words = item.split("-")
            if len(words) != 2:
                raise ValueError("Each item should contain at most one dash.")
            first = int(words[0])
            last = int(words[1])
            if first > last:
                raise ValueError('first=%i > last=%i' % (first, last))
            for number in xrange(first,last+1):
                yield number
        else:
            yield int(item)


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
    '''Iterate over atomic spin multiplicites for the given number of electrons

       **Arguments:**

       nel
            The number of electrons (1-56)

       hund
            When set to True, only one spin multiplicity is returned. Otherwise
            several reasonable spin multiplicities are given.
    '''

    if hund:
        yield mult_presets[nel][0]
    else:
        for mult in mult_presets[nel]:
            yield mult


class Template(BaseTemplate):
    '''A template with modifications to support inclusion of other files.'''
    idpattern = r'[_a-z0-9.:-]+'

    def __init__(self, *args, **kwargs):
        BaseTemplate.__init__(self, *args, **kwargs)
        self.init_include_names()
        self._cache = {}

    def init_include_names(self):
        '''Return a list of include variables

           The include variables in the template are variables of the form
           ${include:name}. This routine lists all the names encountered.
           Duplicates are eliminated.
        '''
        pattern = '%s{(?P<braced>%s)}' % (re.escape(self.delimiter), self.idpattern)
        result = set([])
        for mo in re.finditer(pattern, self.template):
            braced = mo.group('braced')
            if braced is not None and braced.startswith('include:'):
                result.add(braced[8:])
        self.include_names = list(result)

        # log the include names
        if len(self.include_names) > 0 and log.do_medium:
            log('The following includes were detected in the template:')
            for name in self.include_names:
                log('   ', name)

    def load_includes(self, number):
        '''Load included files for a given element number'''
        result = self._cache.setdefault(number, {})
        if len(result) == 0:
            for name in self.include_names:
                with open('%s.%03i' % (name, number)) as f:
                    s = f.read()
                # chop of one final newline if present (mostly the case)
                if s[-1] == '\n':
                    s = s[:-1]
                result['include:%s' % name] = s
        return result


def add_input_arguments(parser):
    '''Add common arguments for any atomdb script.'''
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


epilog_input = '''\
The following fields must be present in the template file: ${charge}, ${mult}
and ${element} or ${number}. It may also contain fields like
${include:some.file} which will be replaced by the contents of the file
"some.file.ZZZ" where ZZZ is the element number left-padded with zeros to fix
the the length at three characters. For example, for oxygen this would be 008.
Such includes may be useful for custom basis sets. The filename of the include
files may only contain characters from the set [_a-z0-9.-].
'''

def add_convert_arguments(parser):
    parser.add_argument('--include-spurious', default=False, action='store_true',
        help='Also convert atoms that are bound by the basis set. These '
             'situations are normally detected by an increase in energy as an '
             'electron is added.')


def iter_states(elements, max_kation, max_anion, hund):
    '''Iterate over all requested atomic states

       **Arguments:**

       elements
            A string that is suitable for ``iter_elements``

       template
            An instance of the ``atomdb.Template`` class

       max_kation
            The limit for the most positive kation

       max_anion
            The limit for the most negative anion

       hund
            Flag to adhere to hund's rule for the spin multiplicities.
    '''
    for number in iter_elements(elements):
        # Loop over all charge states for this element
        for charge in xrange(-max_anion, max_kation+1):
            nel = number - charge
            if nel <= 0:
                continue
            # loop over multiplicities
            for mult in iter_mults(nel, hund):
                yield number, charge, mult


def write_input(number, charge, mult, template, do_overwrite):
    # Directory stuff
    nel = number - charge
    dn_mult = '%03i_%s_%03i_q%+03i/mult%02i' % (
        number, periodic[number].symbol.lower().rjust(2, '_'), nel, charge, mult)
    if not os.path.isdir(dn_mult):
        os.makedirs(dn_mult)

    # Figure out if we want to write
    fn_inp = '%s/atom.in' % dn_mult
    exists = os.path.isfile(fn_inp)
    do_write = not exists or do_overwrite

    if do_write:
        includes = template.load_includes(number)
        with open(fn_inp, 'w') as f:
            f.write(template.substitute(
                includes,
                charge=str(charge),
                mult=str(mult),
                number=str(number),
                element=periodic[number].symbol,
            ))
        if log.do_medium:
            if exists:
                log('Overwritten:      ', fn_inp)
            else:
                log('Written new:      ', fn_inp)
    elif log.do_medium:
        log('Not overwriting:  ', fn_inp)

    return dn_mult, do_write


class EnergyTable(object):
    def __init__(self):
        self.all = {}

    def add(self, number, pop, energy):
        cases = self.all.setdefault(number, {})
        cases[pop] = energy

    def is_spurious(self, number, pop, energy):
        cases = self.all.get(number)
        if cases is not None:
            energy_prev = cases.get(pop-1)
            if energy is not None and energy_prev < energy:
                return True
        return False

    def log(self):
        log(' Nr Pop Chg              Energy          Ionization            Affinity')
        log.hline()
        for number, cases in sorted(self.all.iteritems()):
            for pop, energy in sorted(cases.iteritems()):
                energy_prev = cases.get(pop-1)
                if energy_prev is None:
                    ip_str = ''
                else:
                    ip_str = '% 18.10f' % (energy_prev - energy)

                energy_next = cases.get(pop+1)
                if energy_next is None:
                    ea_str = ''
                else:
                    ea_str = '% 18.10f' % (energy - energy_next)

                log('%3i %3i %+3i  % 18.10f  %18s  %18s' % (
                    number, pop, number-pop, energy, ip_str, ea_str
                ))
            log.blank()
