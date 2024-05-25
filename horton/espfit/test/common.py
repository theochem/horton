# -*- coding: utf-8 -*-
# HORTON: Helpful Open-source Research TOol for N-fermion systems.
# Copyright (C) 2011-2017 The HORTON Development Team
#
# This file is part of HORTON.
#
# HORTON is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
#
# HORTON is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, see <http://www.gnu.org/licenses/>
#
# --


import numpy as np
from horton import Cell
from contextlib import contextmanager


__all__ = ["numpy_seed", "get_random_cell", "check_delta"]


# All, except underflows, is *not* fine.
np.seterr(divide='raise', over='raise', invalid='raise')


@contextmanager
def numpy_seed(seed=1):
    """Temporarily set NumPy's random seed to a given number.

    Parameters
    ----------
    seed : int
           The seed for NumPy's random number generator.
    """
    state = np.random.get_state()
    np.random.seed(seed)
    yield None
    np.random.set_state(state)


def get_random_cell(a, nvec):
    if nvec == 0:
        return Cell(None)
    while True:
        if a <= 0:
            raise ValueError('The first argument must be strictly positive.')
        rvecs = np.random.uniform(0, a, (nvec, 3))
        cell = Cell(rvecs)
        if cell.volume > a**nvec * 0.1:
            return cell


def check_delta(fun, fun_deriv, x, dxs):
    """Check the difference between two function values using the analytical gradient

       Arguments:

       fun
            The function whose derivatives must be to be tested

       fun_deriv
            The implementation of the analytical derivatives

       x
            The argument for the reference point.

       dxs
            A list with small relative changes to x

       For every displacement in ``dxs``, the following computation is repeated:

       1) ``D1 = fun(x+dx) - fun(x)`` is computed.
       2) ``D2 = 0.5*dot(fun_deriv(x+dx) + fun_deriv(x), dx)`` is computed.

       A threshold is set to the median of the D1 set. For each case where |D1|
       is larger than the threshold, |D1 - D2|, should be smaller than the
       threshold.

       This test makes two assumptions:

       1) The gradient at ``x`` is non-zero.
       2) The displacements, ``dxs``, are small enough such that in the majority
          of cases, the linear term in ``fun(x+dx) - fun(x)`` dominates. Hence,
          sufficient elements in ``dxs`` should be provided for this test to
          work.
    """
    assert len(x.shape) == 1
    if len(dxs) < 20:
        raise ValueError('At least 20 displacements are needed for good statistics.')

    dn1s = []
    dn2s = []
    dnds = []
    f0 = fun(x)
    grad0 = fun_deriv(x)
    for dx in dxs:
        f1 = fun(x + dx)
        grad1 = fun_deriv(x + dx)
        grad = 0.5 * (grad0 + grad1)
        d1 = f1 - f0
        if hasattr(d1, '__iter__'):
            norm = np.linalg.norm
        else:
            norm = abs
        d2 = np.dot(grad, dx)

        dn1s.append(norm(d1))
        dn2s.append(norm(d2))
        dnds.append(norm(d1 - d2))
    dn1s = np.array(dn1s)
    dn2s = np.array(dn2s)
    dnds = np.array(dnds)

    # Get the threshold (and mask)
    threshold = np.median(dn1s)
    mask = dn1s > threshold
    # Make sure that all cases for which dn1 is above the treshold, dnd is below
    # the threshold
    if not (dnds[mask] < threshold).all():
        raise AssertionError((
            'The first order approximation on the difference is too wrong. The '
            'threshold is %.1e.\n\nDifferences:\n%s\n\nFirst order '
            'approximation to differences:\n%s\n\nAbsolute errors:\n%s')
            % (threshold,
               ' '.join('%.1e' % v for v in dn1s[mask]),
               ' '.join('%.1e' % v for v in dn2s[mask]),
               ' '.join('%.1e' % v for v in dnds[mask])
               ))
