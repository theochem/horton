#!/usr/bin/env python

from horton import *

# Load the Gaussian output from file from HORTON's test data directory.
fn_fchk = context.get_fn('test/water_sto3g_hf_g03.fchk')
# Replace the previous line with any other fchk file, e.g. fn_fchk = 'yourfile.fchk'.
mol = IOData.from_file(fn_fchk)

# Specify the integration grid
grid = BeckeMolGrid(mol.coordinates, mol.numbers, mol.pseudo_numbers)

# Get the spin-summed density matrix
dm_full = mol.get_dm_full()

# Compute the density on the grid
rho = mol.obasis.compute_grid_density_dm(dm_full, grid.points)

# Compute the expectation value of |r|.
r = (grid.points[:,0]**2 + grid.points[:,1]**2 + grid.points[:,2]**2)**0.5
print 'EXPECTATION VALUE OF |R|:', grid.integrate(rho, r)
