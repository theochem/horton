#!/usr/bin/env python

from horton import *  # pylint: disable=wildcard-import,unused-wildcard-import


# Load the integrals from the file
data = IOData.from_file('hamiltonian_ao.FCIDUMP')

# Access some attributes. In more realistic cases, some code follows that does a
# useful calculation.
print(data.core_energy)
print(data.one_mo[0, 0])


# CODE BELOW IS FOR horton-regression-test.py ONLY. IT IS NOT PART OF THE EXAMPLE.
rt_results = {
    'one_mo': data.one_mo.ravel()[::100],
    'two_mo': data.two_mo.ravel()[::10000],
    'core_energy': data.core_energy,
}
# BEGIN AUTOGENERATED CODE. DO NOT CHANGE MANUALLY.
import numpy as np  # pylint: disable=wrong-import-position
rt_previous = {
    'core_energy': 20.0,
    'one_mo': np.array([
        -98.143299148359503, 0.0, -6.8232742572718648, -0.0014333406428630481,
        -0.0015606243552076162, 0.12476697580764212, 0.0, 0.0
    ]),
    'two_mo': np.array([
        5.9688416466791034, -2.9422895177057246e-07, 0.0, 0.0, -0.00059940882213425917,
        0.0, -1.8851102379070744e-07, 0.0, 0.0, 0.0, 0.0, 1.2785965319495516e-06, 0.0,
        0.0, 1.0046676893268627e-06, 0.0, 0.0, 0.0, 8.6908816309957149e-07, 0.0,
        -2.7850857129776117e-12, 0.0, 0.0, 0.0, 4.8282318941349867e-14, 0.0, 0.0,
        -6.9135313768308591e-07, -2.9422895177057246e-07, 0.0, 0.0, 0.0, 0.0, 0.0,
        -5.6706841716128415e-14, 2.6118860993636095e-24, 1.7350385251530308e-05, 0.0, 0.0,
        0.0, 0.0, 0.0, -3.3218351285338633e-05, 0.0, 0.00012749115201133096, 0.0,
        0.037156902375664713, 0.0, 0.0, -2.3670600726458522e-07, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    ]),
}
