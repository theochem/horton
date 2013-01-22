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
    def __init__(self, number, symbol, cov_radius):
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
            result = self._lookup.get(index.lower())
            if result is None and index.isdigit():
                result = self._lookup.get(int(index))
                if result is None:
                    raise KeyError('Could not find element %s.' % index)
        return result


periodic = Periodic([
    Element(1,'H',0.31*angstrom),
    Element(2,'He',0.28*angstrom),
    Element(3,'Li',1.28*angstrom),
    Element(4,'Be',0.96*angstrom),
    Element(5,'B',0.84*angstrom),
    Element(6,'C',0.73*angstrom),
    Element(7,'N',0.71*angstrom),
    Element(8,'O',0.66*angstrom),
    Element(9,'F',0.57*angstrom),
    Element(10,'Ne',0.58*angstrom),
    Element(11,'Na',1.66*angstrom),
    Element(12,'Mg',1.41*angstrom),
    Element(13,'Al',1.21*angstrom),
    Element(14,'Si',1.11*angstrom),
    Element(15,'P',1.07*angstrom),
    Element(16,'S',1.05*angstrom),
    Element(17,'Cl',1.02*angstrom),
    Element(18,'Ar',1.06*angstrom),
    Element(19,'K',2.03*angstrom),
    Element(20,'Ca',1.76*angstrom),
    Element(21,'Sc',1.70*angstrom),
    Element(22,'Ti',1.60*angstrom),
    Element(23,'V',1.53*angstrom),
    Element(24,'Cr',1.39*angstrom),
    Element(25,'Mn',1.39*angstrom),
    Element(26,'Fe',1.32*angstrom),
    Element(27,'Co',1.26*angstrom),
    Element(28,'Ni',1.24*angstrom),
    Element(29,'Cu',1.32*angstrom),
    Element(30,'Zn',1.22*angstrom),
    Element(31,'Ga',1.22*angstrom),
    Element(32,'Ge',1.20*angstrom),
    Element(33,'As',1.19*angstrom),
    Element(34,'Se',1.20*angstrom),
    Element(35,'Br',1.20*angstrom),
    Element(36,'Kr',1.16*angstrom),
    Element(37,'Rb',2.20*angstrom),
    Element(38,'Sr',1.95*angstrom),
    Element(39,'Y',1.90*angstrom),
    Element(40,'Zr',1.75*angstrom),
    Element(41,'Nb',1.64*angstrom),
    Element(42,'Mo',1.54*angstrom),
    Element(43,'Tc',1.47*angstrom),
    Element(44,'Ru',1.46*angstrom),
    Element(45,'Rh',1.42*angstrom),
    Element(46,'Pd',1.39*angstrom),
    Element(47,'Ag',1.45*angstrom),
    Element(48,'Cd',1.44*angstrom),
    Element(49,'In',1.42*angstrom),
    Element(50,'Sn',1.39*angstrom),
    Element(51,'Sb',1.39*angstrom),
    Element(52,'Te',1.38*angstrom),
    Element(53,'I',1.39*angstrom),
    Element(54,'Xe',1.40*angstrom),
    Element(55,'Cs',2.44*angstrom),
    Element(56,'Ba',2.15*angstrom),
    Element(57,'La',2.07*angstrom),
    Element(58,'Ce',2.04*angstrom),
    Element(59,'Pr',2.03*angstrom),
    Element(60,'Nd',2.01*angstrom),
    Element(61,'Pm',1.99*angstrom),
    Element(62,'Sm',1.98*angstrom),
    Element(63,'Eu',1.98*angstrom),
    Element(64,'Gd',1.96*angstrom),
    Element(65,'Tb',1.94*angstrom),
    Element(66,'Dy',1.92*angstrom),
    Element(67,'Ho',1.92*angstrom),
    Element(68,'Er',1.89*angstrom),
    Element(69,'Tm',1.90*angstrom),
    Element(70,'Yb',1.87*angstrom),
    Element(71,'Lu',1.87*angstrom),
    Element(72,'Hf',1.75*angstrom),
    Element(73,'Ta',1.70*angstrom),
    Element(74,'W',1.62*angstrom),
    Element(75,'Re',1.51*angstrom),
    Element(76,'Os',1.44*angstrom),
    Element(77,'Ir',1.41*angstrom),
    Element(78,'Pt',1.36*angstrom),
    Element(79,'Au',1.36*angstrom),
    Element(80,'Hg',1.32*angstrom),
    Element(81,'Tl',1.45*angstrom),
    Element(82,'Pb',1.46*angstrom),
    Element(83,'Bi',1.48*angstrom),
    Element(84,'Po',1.40*angstrom),
    Element(85,'At',1.50*angstrom),
    Element(86,'Rn',1.50*angstrom),
    Element(87,'Fr',2.60*angstrom),
    Element(88,'Ra',2.21*angstrom),
    Element(89,'Ac',2.15*angstrom),
    Element(90,'Th',2.06*angstrom),
    Element(91,'Pa',2.00*angstrom),
    Element(92,'U',1.96*angstrom),
    Element(93,'Np',1.90*angstrom),
    Element(94,'Pu',1.87*angstrom),
    Element(95,'Am',1.80*angstrom),
    Element(96,'Cm',1.69*angstrom),
    Element(97,'Bk',2.00*angstrom),
    Element(98,'Cf',2.00*angstrom),
    Element(99,'Es',2.00*angstrom),
    Element(100,'Fm',2.00*angstrom),
    Element(101,'Mv',2.00*angstrom),
    Element(102,'No',2.00*angstrom),
    Element(103,'Lr',2.00*angstrom),
    Element(104,'Rf',2.00*angstrom),
    Element(105,'Db',2.00*angstrom),
    Element(106,'Sg',2.00*angstrom),
    Element(107,'Bh',2.00*angstrom),
    Element(108,'Hs',2.00*angstrom),
    Element(109,'Mt',2.00*angstrom),
    Element(110,'Uun',2.00*angstrom),
    Element(111,'Uuu',2.00*angstrom),
    Element(112,'Uub',2.00*angstrom),
    Element(113,'Uut',2.00*angstrom),
    Element(114,'Uuq',2.00*angstrom),
    Element(115,'Uup',2.00*angstrom),
    Element(116,'Uuh',2.00*angstrom),
    Element(117,'Uus',2.00*angstrom),
    Element(118,'Uuo',2.00*angstrom),
])
