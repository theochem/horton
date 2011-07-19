# -*- coding: utf-8 -*-
# Horton is a Density Functional Theory program.
# Copyright (C) 2011 Toon Verstraelen <Toon.Verstraelen@UGent.be>, ...
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


class Element(object):
    '''Represents an element from the periodic table.

       The following attributes are supported:

       number
            The atomic number

       symbol
            A string with the symbol of the element.
    '''
    def __init__(self, number, symbol):
        self.number = number
        self.symbol = symbol


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
    Element(1,'H'),
    Element(2,'He'),
    Element(3,'Li'),
    Element(4,'Be'),
    Element(5,'B'),
    Element(6,'C'),
    Element(7,'N'),
    Element(8,'O'),
    Element(9,'F'),
    Element(10,'Ne'),
    Element(11,'Na'),
    Element(12,'Mg'),
    Element(13,'Al'),
    Element(14,'Si'),
    Element(15,'P'),
    Element(16,'S'),
    Element(17,'Cl'),
    Element(18,'Ar'),
    Element(19,'K'),
    Element(20,'Ca'),
    Element(21,'Sc'),
    Element(22,'Ti'),
    Element(23,'V'),
    Element(24,'Cr'),
    Element(25,'Mn'),
    Element(26,'Fe'),
    Element(27,'Co'),
    Element(28,'Ni'),
    Element(29,'Cu'),
    Element(30,'Zn'),
    Element(31,'Ga'),
    Element(32,'Ge'),
    Element(33,'As'),
    Element(34,'Se'),
    Element(35,'Br'),
    Element(36,'Kr'),
    Element(37,'Rb'),
    Element(38,'Sr'),
    Element(39,'Y'),
    Element(40,'Zr'),
    Element(41,'Nb'),
    Element(42,'Mo'),
    Element(43,'Tc'),
    Element(44,'Ru'),
    Element(45,'Rh'),
    Element(46,'Pd'),
    Element(47,'Ag'),
    Element(48,'Cd'),
    Element(49,'In'),
    Element(50,'Sn'),
    Element(51,'Sb'),
    Element(52,'Te'),
    Element(53,'I'),
    Element(54,'Xe'),
    Element(55,'Cs'),
    Element(56,'Ba'),
    Element(57,'La'),
    Element(58,'Ce'),
    Element(59,'Pr'),
    Element(60,'Nd'),
    Element(61,'Pm'),
    Element(62,'Sm'),
    Element(63,'Eu'),
    Element(64,'Gd'),
    Element(65,'Tb'),
    Element(66,'Dy'),
    Element(67,'Ho'),
    Element(68,'Er'),
    Element(69,'Tm'),
    Element(70,'Yb'),
    Element(71,'Lu'),
    Element(72,'Hf'),
    Element(73,'Ta'),
    Element(74,'W'),
    Element(75,'Re'),
    Element(76,'Os'),
    Element(77,'Ir'),
    Element(78,'Pt'),
    Element(79,'Au'),
    Element(80,'Hg'),
    Element(81,'Tl'),
    Element(82,'Pb'),
    Element(83,'Bi'),
    Element(84,'Po'),
    Element(85,'At'),
    Element(86,'Rn'),
    Element(87,'Fr'),
    Element(88,'Ra'),
    Element(89,'Ac'),
    Element(90,'Th'),
    Element(91,'Pa'),
    Element(92,'U'),
    Element(93,'Np'),
    Element(94,'Pu'),
    Element(95,'Am'),
    Element(96,'Cm'),
    Element(97,'Bk'),
    Element(98,'Cf'),
    Element(99,'Es'),
    Element(100,'Fm'),
    Element(101,'Mv'),
    Element(102,'No'),
    Element(103,'Lr'),
    Element(104,'Rf'),
    Element(105,'Db'),
    Element(106,'Sg'),
    Element(107,'Bh'),
    Element(108,'Hs'),
    Element(109,'Mt'),
    Element(110,'Uun'),
    Element(111,'Uuu'),
    Element(112,'Uub'),
    Element(113,'Uut'),
    Element(114,'Uuq'),
    Element(115,'Uup'),
    Element(116,'Uuh'),
    Element(117,'Uus'),
    Element(118,'Uuo'),
])
