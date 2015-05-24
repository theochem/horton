#!/usr/bin/env python

from horton import *

# Define number of occupied orbitals and total number of basis functions
nocc = 23
nbasis = 115

# Define linalg factory
lf = DenseLinalgFactory(nbasis)

# Read Hamiltonian from file 'FCIDUMP'
one, two, core = load_fcidump(lf, 'FCIDUMP')
