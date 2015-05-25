#!/usr/bin/env python

from horton import *

# Load the integrals from the file
mol = Molecule.from_file('hamiltonian_ao.h5')

# Access some attributes. In more realistic cases, some code follows that does a
# useful calculation.
print mol.kin.nbasis
