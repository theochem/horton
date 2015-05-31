#!/usr/bin/env python

import numpy as np
from horton import *

# define the ring
natom = 11
spacing = 1.3 # distance between two neighboring atoms in bohr
radius = spacing/(2*np.sin(np.pi/natom))

# define the coordinates and elements
coordinates = np.zeros((natom, 3))
numbers = np.ones(natom, dtype=int) # must be integers
for iatom in xrange(natom):
    angle = (2*np.pi/natom)*iatom
    coordinates[iatom, 0] = radius*np.cos(angle)
    coordinates[iatom, 1] = radius*np.sin(angle)

# write the molecule to an XYZ file (optional)
mol = IOData(coordinates=coordinates, numbers=numbers, title='H Ring')
mol.to_file('ring.xyz')
