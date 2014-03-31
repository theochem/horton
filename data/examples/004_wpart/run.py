#!/usr/bin/env python

import numpy as np
from horton import *

# Load the Gaussian output from file
fn_fchk = context.get_fn('test/water_sto3g_hf_g03.fchk')
# Replace the previous line with any other fchk file, e.g. fn_fchk = 'yourfile.fchk'.
mol = Molecule.from_file(fn_fchk)

# Partition the density with the Becke scheme
grid = BeckeMolGrid(mol.coordinates, mol.numbers, mol.pseudo_numbers, mode='only')
moldens = mol.obasis.compute_grid_density_dm(mol.wfn.dm_full, grid.points)
bp = BeckeWPart(mol.coordinates, mol.numbers, mol.pseudo_numbers, grid, moldens, local=True)
bp.do_charges()

# Write the result to a file
np.savetxt('charges.txt', bp['charges'])
