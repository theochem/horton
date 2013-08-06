#!/usr/bin/env python

import numpy as np
from horton import *

# Load the Gaussian output from file
fn_fchk = context.get_fn('test/water_sto3g_hf_g03.fchk')
# Replace the previous line with any other fchk file, e.g. fn_fchk = 'yourfile.fchk'.
sys = System.from_file(fn_fchk)

# Partition the density with the Becke scheme
grid = BeckeMolGrid(sys, mode='keep')
bp = BeckeWPart(sys, grid, local=True)
bp.do_charges()

# Write the result to a file
np.savetxt('charges.txt', bp['charges'])
