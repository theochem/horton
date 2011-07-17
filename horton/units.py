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
"""Conversion from and to atomic units

   Internally the Horton package works always in atomic units. This unit system
   is consistent, like the SI unit system one does not need conversion factors
   in the middle of a computation once all values are converted to atomic units.
   This facilitates the programming and reduces accidental bugs due to
   forgetting these conversion factor in the body of the code.

   References for the conversion values:

   * B. J. Mohr and B. N. Taylor,
     CODATA recommended values of the fundamental physical
     constants: 1998, Rev. Mod. Phys. 72(2), 351 (2000)
   * The NIST Reference on Constants, Units, and Uncertainty
     (http://physics.nist.gov/cuu/Constants/index.html)
   * 1 calorie = 4.184 Joules

   Naming conventions in this module: unit is the value of one external unit
   in internal - i.e. atomic - units. e.g. If you want to have a distance of
   five angstrom in internal units: ``5*angstrom``. If you want to convert a
   length of 5 internal units to angstrom: ``5/angstrom``. It is recommended to
   perform this kind of conversions, only when data is read from the input and
   data is written to the output.

   An often recurring question is how to convert a frequency in internal units
   to a spectroscopic wavenumber in inverse centimeters. This is how it can be
   done::

     >>> from horton import centimeter, lightspeed
     >>> invcm = lightspeed/centimeter
     >>> freq = 0.00320232
     >>> print freq/invcm

   These are the conversion constants defined in this module:

"""


from horton.constants import avogadro


# *** Generic ***
au = 1.0


# *** Charge ***

coulomb = 1.0/1.602176462e-19

# Mol

mol = avogadro

# *** Mass ***

kilogram = 1.0/9.10938188e-31

gram = 1.0e-3*kilogram
miligram = 1.0e-6*kilogram
unified = 1.0e-3*kilogram/mol
amu = unified

# *** Length ***

meter = 1.0/0.5291772083e-10

decimeter = 1.0e-1*meter
centimeter = 1.0e-2*meter
milimeter = 1.0e-3*meter
micrometer = 1.0e-6*meter
nanometer = 1.0e-9*meter
angstrom = 1.0e-10*meter
picometer = 1.0e-12*meter

# *** Volume ***

liter = decimeter**3

# *** Energy ***

joule = 1/4.35974381e-18

calorie = 4.184*joule
kjmol = 1.0e3*joule/mol
kcalmol = 1.0e3*calorie/mol
electronvolt = (1.0/coulomb)*joule
rydberg = 0.5

# *** Force ***

newton = joule/meter

# *** Angles ***

deg = 0.017453292519943295
rad = 1.0

# *** Time ***

second = 1/2.418884326500e-17

nanosecond = 1e-9*second
femtosecond = 1e-15*second
picosecond = 1e-12*second

# *** Frequency ***

hertz = 1/second

# *** Pressure ***

pascal = newton/meter**2
bar = 100000*pascal
atm = 1.01325*bar

# *** Temperature ***

kelvin = 1.0

# *** Dipole ***

debye = 0.39343031369146675 # = 1e-21*coulomb*meter**2/second/lightspeed

# *** Current ***

ampere = coulomb/second
