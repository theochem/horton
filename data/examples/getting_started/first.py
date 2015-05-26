#!/usr/bin/env python

from horton import *

# Load the molecule from an ``.xyz`` file.
mol = Molecule.from_file(context.get_fn('test/water.xyz'))

# Compute the molecular mass
mass = 0.0
#   Loop over all atomic numbers
for number in mol.numbers:
    mass += periodic[number].mass

# Print the mass in the amu unit
print 'MOLECULAR MASS [amu]: %.5f' % (mass/amu)

# Store the mass in the Molecule instance, in order to write it to the file
mol.mass = mass

# Write data in the mol object to a file in Horton's internal HDF5-based
# file format.
mol.to_file('water.h5')
