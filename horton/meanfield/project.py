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
"""Projection of 1-electron orbitals to a new basis set"""

# TODO: Move to gobasis?

import numpy as np
from scipy.linalg import sqrtm

from horton.gbasis.cext import GOBasis

__all__ = ['ProjectionError', 'project_orbitals_mgs', 'project_orbitals_ortho']


class ProjectionError(Exception):
    pass


def project_orbitals_mgs(obasis0, obasis1, orb0, orb1, eps=1e-10):
    """Project the orbitals onto a new basis set with the modified Gram-Schmidt algorithm.

    The orbitals in ``orb0`` (w.r.t. ``obasis0``) are projected onto ``obasis1`` and
    stored in ``orb1``.

    Parameters
    ----------

    obasis0 : GOBasis
             The orbital basis for the original wavefunction expansion.
    obasis1 : GOBasis
             The new orbital basis for the projected wavefunction expansion.
    orb0 : Orbitals
          The expansion of the original orbitals.
    orb1 : Orbitals
          An output argument in which the projected orbitals will be stored.
    eps : float
         A threshold for the renormalization in the Gram-Schmidt procedure

    Notes
    -----

    The projection is based on the Modified Gram-Schmidt (MGS) process. In
    each iteration of the MGS, a renormalization is carried out. If the norm
    in this step is smaller than ``eps``, an error is raised.

    Note that ``orb1`` will be incomplete in several ways. The orbital
    energies are not copied. Only the occupied orbitals in ``orb0`` are
    projected. Coefficients of higher orbitals are set to zero. The orbital
    occupations are simply copied. This should be sufficient to construct
    an initial guess in a new orbital basis set based on a previous solution.

    If the number of orbitals in ``orb1`` is too small to store all projected
    orbitals, an error is raised.
    """
    # Compute the overlap matrix of the combined orbital basis
    obasis_both = GOBasis.concatenate(obasis0, obasis1)
    olp_both = obasis_both.compute_overlap()

    # Select the blocks of interest from the big overlap matrix
    olp_21 = olp_both[obasis0.nbasis:, :obasis0.nbasis]
    olp_22 = olp_both[obasis0.nbasis:, obasis0.nbasis:]

    # Construct the projector. This minimizes the L2 norm between the new and old
    # orbitals, which does not account for orthonormality.
    projector = np.dot(np.linalg.pinv(olp_22), olp_21)

    # Project occupied orbitals.
    i1 = 0
    for i0 in xrange(orb0.nfn):
        if orb0.occupations[i0] == 0.0:
            continue
        if i1 > orb1.nfn:
            raise ProjectionError('Not enough functions available in orb1 to store the '
                                  'projected orbitals.')
        orb1.coeffs[:, i1] = np.dot(projector, orb0.coeffs[:, i0])
        orb1.occupations[i1] = orb0.occupations[i0]
        i1 += 1

    # clear all parts of orb1 that were not touched by the projection loop
    ntrans = i1
    del i1
    orb1.coeffs[:, ntrans:] = 0.0
    orb1.occupations[ntrans:] = 0.0
    orb1.energies[:] = 0.0

    # auxiliary function for the MGS algo
    def dot22(a, b):
        return np.dot(np.dot(a, olp_22), b)

    # Apply the MGS algorithm to orthogonalize the orbitals
    for i1 in xrange(ntrans):
        orb = orb1.coeffs[:, i1]

        # Subtract overlap with previous orbitals
        for j1 in xrange(i1):
            other = orb1.coeffs[:, j1]
            orb -= other * dot22(other, orb) / np.sqrt(dot22(other, other))

        # Renormalize
        norm = np.sqrt(dot22(orb, orb))
        if norm < eps:
            raise ProjectionError('The norm of a vector in the MGS algorithm becomes too '
                                  'small. Orbitals are redundant in new basis.')
        orb /= norm


def project_orbitals_ortho(olp0, olp1, orb0, orb1):
    r"""Re-orthogonalize the orbitals .

    The orbitals in ``orb0`` (w.r.t. ``obasis0``) are re-orthonormalized w.r.t.
    ``obasis1`` and stored in ``orb1``.

    Parameters
    ----------

    olp0 : TwoIndex or GOBasis
           The overlap matrix (or alternatively the orbital basis) for the original
           wavefunction expansion.
    olp1 : TwoIndex or GOBasis
           The overlap matrix (or alternatively the orbital basis) for the projected
           wavefunction expansion.
    orb0 : DenseExpansion
           The expansion of the original orbitals.
    orb1 : DenseExpansion
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
    """

    def helper_olp(olp):
        if isinstance(olp, GOBasis):
            obasis = olp
            olp = obasis.compute_overlap()
        elif isinstance(olp, np.ndarray):
            obasis = None
        else:
            raise TypeError('The olp arguments must be an instance of TwoIndex or GOBasis.')
        return olp, obasis

    olp0, obasis0 = helper_olp(olp0)
    olp1, obasis1 = helper_olp(olp1)

    # Transform the coefficients
    tf = sqrtm(np.dot(np.linalg.inv(olp1), olp0))
    orb1.coeffs[:] = np.dot(tf, orb0.coeffs)
    # Clear the energies in orb1 as they can not be defined in a meaningful way
    orb1.energies[:] = 0.0
    # Just copy the occupation numbers and hope for the best
    orb1.occupations[:] = orb0.occupations
