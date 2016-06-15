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
'''Conversion from and to atomic units

   Internally HORTON always uses atomic units. Atomic units are consistent,
   similar to the SI unit system: one does not need conversion factors in the
   middle of a computation. This choice facilitates the programming and reduces
   accidental bugs.

   References for the conversion values:

   * B. J. Mohr and B. N. Taylor,
     CODATA recommended values of the fundamental physical
     constants: 1998, Rev. Mod. Phys. 72(2), 351 (2000)
   * The NIST Reference on Constants, Units, and Uncertainty
     (http://physics.nist.gov/cuu/Constants/index.html)
   * 1 calorie = 4.184 Joules

   **Conventions followed by this module:**

   Let foo be is the value of an external unit in internal (atomic) units. The
   way to use this unit is as follows: ``5*foo`` litterally means `five times
   foo`. The result of this operation is a floating point number for this value
   in atomic units.

   **Examples:**

   If you want to have a distance of five angstrom in internal units:
   ``5*angstrom``.

   If you want to convert a length of 5 internal units to angstrom:
   ``5/angstrom``.

   **Remarks:**

   It is highly recommended to perform unit conversions only when data is read
   from the input or data is written to the output. It may also be useful in
   `input scripts` that use HORTON. Do not perform any unit conversion in other
   parts of the program.

   An often recurring question is how to convert a frequency in internal units
   to a spectroscopic wavenumber in inverse centimeters. This is how it can be
   done::

     >>> from horton import centimeter, lightspeed
     >>> invcm = lightspeed/centimeter
     >>> freq = 0.00320232
     >>> print freq/invcm

   These are the conversion constants defined in this module:

'''


from horton.constants import avogadro


# *** Generic ***
au = 1.0


# *** Charge ***

coulomb = 1.0/1.602176462e-19

# *** Mass ***

kilogram = 1.0/9.10938188e-31

gram = 1.0e-3*kilogram
miligram = 1.0e-6*kilogram
unified = 1.0e-3*kilogram/avogadro
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
kjmol = 1.0e3*joule/avogadro
kcalmol = 1.0e3*calorie/avogadro
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



# automatically spice up the docstrings

lines = [
    '    ================  ==================',
    '    Name              Value             ',
    '    ================  ==================',
]

for key, value in sorted(globals().iteritems()):
    if not isinstance(value, float):
        continue
    lines.append('    %16s  %.10e' % (key, value))
lines.append('    ================  ==================')

__doc__ += '\n'.join(lines)


del lines
