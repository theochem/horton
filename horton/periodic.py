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

from horton.context import *
import numpy as np


'''Periodic table of elements.

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


__all__ = ['periodic', 'Element', 'Periodic']


from horton.units import angstrom


class Element(object):
    '''Represents an element from the periodic table.

       The following attributes are supported:

       number
            The atomic number

       symbol
            A string with the symbol of the element.

       cov_radius
            The covalent radius
    '''
    def __init__(self, number= None, symbol= None, cov_radius=None):
        self.number = number
        self.symbol = symbol
        self.cov_radius = cov_radius


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
        'str': (lambda s: s.strip()),
        'angstrom': (lambda s: float(s)*angstrom),
    }

    filename = 'elements'
    fn = context.get_fn('%s.txt' % filename)
    nelement = 0
    with open(fn,'r') as infile:
        line1= infile.readline()
        number = int(line1.split(':')[1])
        rows = {}
        step = 1
        for line in infile:
            line = line[:line.find('#')].strip()
            if len(line) > 0:
                if step == 1:
                    name = line.split(',')[0]
                    convert = line.split(',')[1]
                    step = 2
                elif step == 2 :
                    row = [convertors[convert](l) for l in line.split(",")]
                    nelement = max(nelement, len(row))
                    rows[name] = row
                    step = 1

    elements=[]
    for i in xrange(nelement):
        kwargs = dict((name, values[i]) for name, values in rows.iteritems())
        elements.append(Element(**kwargs))


    return Periodic(elements)


periodic = load_periodic()
