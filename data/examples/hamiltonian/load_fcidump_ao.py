#!/usr/bin/env python

from horton import *

# Load the integrals from the file
data = IOData.from_file('hamiltonian_ao.FCIDUMP')

# Access some attributes. In more realistic cases, some code follows that does a
# useful calculation.
print data.core_energy
print data.one_mo.get_element(0, 0)

# Assign results to variables for regression testing
# --------------------------------------------------
result_energy = data.core_energy
# --------------------------------------------------
