# -*- coding: utf-8 -*-
# Horton is a development platform for electronic structure methods.
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
'''Periodic table of elements

   This module contains an object ``periodic`` that can be used as a Pythonic
   periodic table. It can be used as follows::

       >>> from horton import periodic
       >>> periodic['si'].number
       14
       >>> periodic['He'].number
       2
       >>> periodic['h'].symbol
       'H'
       >>> periodic[3].symbol
       'Li'
       >>> periodic['5'].symbol
       'B'
'''


from horton.context import context
from horton.units import angstrom


__all__ = ['periodic', 'Element', 'Periodic']



class Element(object):
    '''Represents an element from the periodic table.

       The following attributes are supported:

       number
            The atomic number

       symbol
            A string with the symbol of the element.

       cov_radius
            The covalent radius. B. Cordero, V. Gomez, A. E. Platero-Prats, M.
            Reves, J. Echeverria, E. Cremades, F. Barragan, and S. Alvarez,
            Dalton Trans. pp. 2832--2838 (2008), URL
            http://dx.doi.org/10.1039/b801115j

       bs_radius
            The Bragg-Slater radius. J. C. Slater, J. Chem. Phys. 41, 3199
            (1964), URL http://dx.doi.org/10.1063/1.1725697

       vdw_radius
            van der Waals radius. R. S. Rowland and R. Taylor, J. Phys. Chem.
            100, 7384 (1996), URL http://dx.doi.org/10.1021/jp953141+

       wc_radius
            Waber-Cromer radius of the outermost orbital maximum. J. T. Waber
            and D. T. Cromer, J. Chem. Phys. 42, 4116 (1965), URL
            http://dx.doi.org/10.1063/1.1695904
    '''
    def __init__(self, number=None, symbol=None, cov_radius=None, bs_radius=None, vdw_radius=None, wc_radius=None):
        self.number = number
        self.symbol = symbol
        self.cov_radius = cov_radius
        self.bs_radius = bs_radius
        self.vdw_radius = vdw_radius
        self.wc_radius = wc_radius


class Periodic(object):
    '''A periodic table data structure.'''
    def __init__(self, elements):
        '''**Arguments:**

           elements
                A list of :class:`Element` instances.
        '''
        self.elements = elements
        self._lookup = {}
        for element in elements:
            self._lookup[element.number] = element
            self._lookup[element.symbol.lower()] = element

    def __getitem__(self, index):
        '''Get an element from the table based on a flexible index.

           **Argument:**

           index
                This can be either an integer atomic number, a string with the
                elemental symbol (any case), or a string with the atomic number.

           **Returns:** the corresponding :class:`Element` instance
        '''
        result = self._lookup.get(index)
        if result is None and isinstance(index, basestring):
            index = index.strip()
            result = self._lookup.get(index.lower())
            if result is None and index.isdigit():
                result = self._lookup.get(int(index))
                if result is None:
                    raise KeyError('Could not find element %s.' % index)
        return result

def load_periodic():

    convertors = {
        'int': (lambda s: int(s)),
        'float': (lambda s : float(s)),
        'str': (lambda s: s.strip()),
        'angstrom': (lambda s: float(s)*angstrom),
    }

    fn = context.get_fn('elements.txt')
    nelement = 0
    with open(fn,'r') as infile:
        rows = {}
        step = 1
        for line in infile:
            line = line[:line.find('#')].strip()
            if len(line) > 0:
                if step == 1:
                    name, convert = line.split()
                    step = 2
                elif step == 2:
                    row = []
                    for word in line.split():
                        if word == 'None':
                            row.append(None)
                        else:
                            row.append(convertors[convert](word))
                    nelement = max(nelement, len(row))
                    rows[name] = row
                    step = 1

    elements=[]
    args=[]
    for i in xrange(nelement):
        for name, values in rows.iteritems():
            if len(values) < nelement :
                for j in range(nelement-len(values)):
                    values.append(None)
            args.append((name,values[i]))
        kwargs = dict(args)
        elements.append(Element(**kwargs))


    return Periodic(elements)


periodic = load_periodic()
