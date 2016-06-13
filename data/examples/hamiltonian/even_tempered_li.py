#!/usr/bin/env python

import numpy as np
from horton import *

# specify the even tempered basis set
alpha_low = 5e-3
alpha_high = 5e2
nbasis = 30
lnratio = (np.log(alpha_high) - np.log(alpha_low))/(nbasis-1)

# build a list of "contractions". These aren't real contractions as every
# contraction only contains one basis function.
bcs = []
for ibasis in xrange(nbasis):
    alpha = alpha_low * np.exp(lnratio * ibasis)
    # arguments of GOBasisContraction:
    #     shell_type, list of exponents, list of contraction coefficients
    bcs.append(GOBasisContraction(0, np.array([alpha]), np.array([1.0])))

# Finish setting up the basis set:
ba = GOBasisAtom(bcs)
obasis = get_gobasis(np.array([[0.0, 0.0, 0.0]]), np.array([3]), default=ba)
