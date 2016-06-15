# -*- coding: utf-8 -*-
# HORTON: Helpful Open-source Research TOol for N-fermion systems.
# Copyright (C) 2011-2016 The HORTON Development Team
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
'''Projection of 1-electron orbitals to a new basis set'''


import numpy as np
from horton.matrix.base import TwoIndex
from horton.matrix.dense import DenseLinalgFactory
from horton.gbasis.cext import GOBasis


__all__ = ['ProjectionError', 'project_orbitals_mgs', 'project_orbitals_ortho']


class ProjectionError(Exception):
    pass


def project_orbitals_mgs(obasis0, obasis1, exp0, exp1, eps=1e-10):
    '''Project the orbitals onto a new basis set with the modified Gram-Schmidt algorithm.

    The orbitals in ``exp0`` (w.r.t. ``obasis0``) are projected onto ``obasis1`` and
    stored in ``exp1``.

    Parameters
    ----------

    obasis0 : GOBasis
             The orbital basis for the original wavefunction expansion.
    obasis1 : GOBasis
             The new orbital basis for the projected wavefunction expansion.
    exp0 : DenseExpansion
          The expansion of the original orbitals.
    exp1 : DenseExpansion
          An output argument in which the projected orbitals will be stored.
    eps : float
         A threshold for the renormalization in the Gram-Schmidt procedure

    Notes
    -----

    The projection is based on the Modified Gram-Schmidt (MGS) process. In
    each iteration of the MGS, a renormalization is carried out. If the norm
    in this step is smaller than ``eps``, an error is raised.

    Note that ``exp1`` will be incomplete in several ways. The orbital
    energies are not copied. Only the occupied orbitals in ``exp0`` are
    projected. Coefficients of higher orbitals are set to zero. The orbital
    occupations are simply copied. This should be sufficient to construct
    an initial guess in a new orbital basis set based on a previous solution.

    If the number of orbitals in ``exp1`` is too small to store all projected
    orbitals, an error is raised.
    '''
    # Compute the overlap matrix of the combined orbital basis
    obasis_both = GOBasis.concatenate(obasis0, obasis1)
    lf = DenseLinalgFactory(obasis_both.nbasis)
    olp_both = obasis_both.compute_overlap(lf)

    # Select the blocks of interest from the big overlap matrix
    olp_21 = olp_both._array[obasis0.nbasis:, :obasis0.nbasis]
    olp_22 = olp_both._array[obasis0.nbasis:, obasis0.nbasis:]

    # Construct the projector. This minimizes the L2 norm between the new and old
    # orbitals, which does not account for orthonormality.
    projector = np.dot(np.linalg.pinv(olp_22), olp_21)

    # Project occupied orbitals.
    i1 = 0
    for i0 in xrange(exp0.nfn):
        if exp0.occupations[i0] == 0.0:
            continue
        if i1 > exp1.nfn:
            raise ProjectionError('Not enough functions available in exp1 to store the '
                                  'projected orbitals.')
        exp1.coeffs[:,i1] = np.dot(projector, exp0.coeffs[:,i0])
        exp1.occupations[i1] = exp0.occupations[i0]
        i1 += 1

    # clear all parts of exp1 that were not touched by the projection loop
    ntrans = i1
    del i1
    exp1.coeffs[:,ntrans:] = 0.0
    exp1.occupations[ntrans:] = 0.0
    exp1.energies[:] = 0.0

    # auxiliary function for the MGS algo
    def dot22(a, b):
        return np.dot(np.dot(a, olp_22), b)

    # Apply the MGS algorithm to orthogonalize the orbitals
    for i1 in xrange(ntrans):
        orb = exp1.coeffs[:, i1]

        # Subtract overlap with previous orbitals
        for j1 in xrange(i1):
            other = exp1.coeffs[:, j1]
            orb -= other*dot22(other, orb)/np.sqrt(dot22(other, other))

        # Renormalize
        norm = np.sqrt(dot22(orb, orb))
        if norm < eps:
            raise ProjectionError('The norm of a vector in the MGS algorithm becomes too '
                                  'small. Orbitals are redundant in new basis.')
        orb /= norm


def project_orbitals_ortho(olp0, olp1, exp0, exp1):
    r'''Re-orthogonalize the orbitals .

    The orbitals in ``exp0`` (w.r.t. ``obasis0``) are re-orthonormalized w.r.t.
    ``obasis1`` and stored in ``exp1``.

    Parameters
    ----------

    olp0 : TwoIndex or GOBasis
           The overlap matrix (or alternatively the orbital basis) for the original
           wavefunction expansion.
    olp1 : TwoIndex or GOBasis
           The overlap matrix (or alternatively the orbital basis) for the projected
           wavefunction expansion.
    exp0 : DenseExpansion
           The expansion of the original orbitals.
    exp1 : DenseExpansion
           An output argument in which the projected orbitals will be stored.

    Notes
    -----

    This projection just transforms the old orbitals to an orthogonal basis
    by a multiplication with the square root of the old overlap matrix. The
    orbitals in this basis are then again multiplied with the inverse square
    root of the new overlap matrix:

    .. math ::

        C_\text{new} = S_\text{new}^{-1/2} S_\text{old}^{1/2} C_\text{old}

    This guarantees that :math:`C_\text{new}^T S_\text{new} C_\text{new} = I`
    if :math:`C_\text{old}^T S_\text{old} C_\text{old} = I`. This approach is
    simple and robust but the current implementation has some limitations: it
    only works for projections between basis sets of the same size and it
    assumes that there is some similarity between the new and old
    orthogonalized atomic basis sets. The latter is only the case when the
    old and new atomic basis sets are very similar, e.g. for small geometric
    changes.
    '''
    if olp0.nbasis != olp1.nbasis:
        raise ValueError('The two basis sets must have the same size')

    def helper_olp(olp):
        if isinstance(olp, GOBasis):
            obasis = olp
            lf = DenseLinalgFactory(obasis.nbasis)
            olp = obasis.compute_overlap(lf)
        elif isinstance(olp, TwoIndex):
            obasis = None
        else:
            raise TypeError('The olp arguments must be an instance of TwoIndex or GOBasis.')
        return olp, obasis

    olp0, obasis0 = helper_olp(olp0)
    olp1, obasis1 = helper_olp(olp1)

    tmp = olp1.inverse()
    tmp.idot(olp0)
    tf = tmp.sqrt()

    # Transform the coefficients
    exp1.assign_dot(tf, exp0)

    # Clear the energies in exp1 as they can not be defined in a meaningful way
    exp1.energies[:] = 0.0
    # Just copy the occupation numbers and hope for the best
    exp1.occupations[:] = exp0.occupations
