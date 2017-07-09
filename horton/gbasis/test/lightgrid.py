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


__all__ = ['generate_molecular_grid', 'integrate']


lebedev_50_points = np.array([
    [1.0, 0.0, 0.0],
    [-1.0, 0.0, 0.0],
    [0.0, 1.0, 0.0],
    [0.0, -1.0, 0.0],
    [0.0, 0.0, 1.0],
    [0.0, 0.0, -1.0],
    [0.0, 0.7071067811865475, 0.7071067811865475],
    [0.0, 0.7071067811865475, -0.7071067811865475],
    [0.0, -0.7071067811865475, 0.7071067811865475],
    [0.0, -0.7071067811865475, -0.7071067811865475],
    [0.7071067811865475, 0.0, 0.7071067811865475],
    [0.7071067811865475, 0.0, -0.7071067811865475],
    [-0.7071067811865475, 0.0, 0.7071067811865475],
    [-0.7071067811865475, 0.0, -0.7071067811865475],
    [0.7071067811865475, 0.7071067811865475, 0.0],
    [0.7071067811865475, -0.7071067811865475, 0.0],
    [-0.7071067811865475, 0.7071067811865475, 0.0],
    [-0.7071067811865475, -0.7071067811865475, 0.0],
    [0.5773502691896258, 0.5773502691896258, 0.5773502691896258],
    [0.5773502691896258, 0.5773502691896258, -0.5773502691896258],
    [0.5773502691896258, -0.5773502691896258, 0.5773502691896258],
    [0.5773502691896258, -0.5773502691896258, -0.5773502691896258],
    [-0.5773502691896258, 0.5773502691896258, 0.5773502691896258],
    [-0.5773502691896258, 0.5773502691896258, -0.5773502691896258],
    [-0.5773502691896258, -0.5773502691896258, 0.5773502691896258],
    [-0.5773502691896258, -0.5773502691896258, -0.5773502691896258],
    [0.3015113445777636, 0.3015113445777636, 0.9045340337332909],
    [0.3015113445777636, 0.3015113445777636, -0.9045340337332909],
    [0.3015113445777636, -0.3015113445777636, 0.9045340337332909],
    [0.3015113445777636, -0.3015113445777636, -0.9045340337332909],
    [-0.3015113445777636, 0.3015113445777636, 0.9045340337332909],
    [-0.3015113445777636, 0.3015113445777636, -0.9045340337332909],
    [-0.3015113445777636, -0.3015113445777636, 0.9045340337332909],
    [-0.3015113445777636, -0.3015113445777636, -0.9045340337332909],
    [0.3015113445777636, 0.9045340337332909, 0.3015113445777636],
    [0.3015113445777636, -0.9045340337332909, 0.3015113445777636],
    [0.3015113445777636, 0.9045340337332909, -0.3015113445777636],
    [0.3015113445777636, -0.9045340337332909, -0.3015113445777636],
    [-0.3015113445777636, 0.9045340337332909, 0.3015113445777636],
    [-0.3015113445777636, -0.9045340337332909, 0.3015113445777636],
    [-0.3015113445777636, 0.9045340337332909, -0.3015113445777636],
    [-0.3015113445777636, -0.9045340337332909, -0.3015113445777636],
    [0.9045340337332909, 0.3015113445777636, 0.3015113445777636],
    [-0.9045340337332909, 0.3015113445777636, 0.3015113445777636],
    [0.9045340337332909, 0.3015113445777636, -0.3015113445777636],
    [-0.9045340337332909, 0.3015113445777636, -0.3015113445777636],
    [0.9045340337332909, -0.3015113445777636, 0.3015113445777636],
    [-0.9045340337332909, -0.3015113445777636, 0.3015113445777636],
    [0.9045340337332909, -0.3015113445777636, -0.3015113445777636],
    [-0.9045340337332909, -0.3015113445777636, -0.3015113445777636]])

lebedev_50_weights = np.array([
    0.0126984126984127, 0.0126984126984127, 0.0126984126984127, 0.0126984126984127,
    0.0126984126984127, 0.0126984126984127, 0.02257495590828924, 0.02257495590828924,
    0.02257495590828924, 0.02257495590828924, 0.02257495590828924, 0.02257495590828924,
    0.02257495590828924, 0.02257495590828924, 0.02257495590828924, 0.02257495590828924,
    0.02257495590828924, 0.02257495590828924, 0.02109375, 0.02109375, 0.02109375,
    0.02109375, 0.02109375, 0.02109375, 0.02109375, 0.02109375, 0.02017333553791887,
    0.02017333553791887, 0.02017333553791887, 0.02017333553791887, 0.02017333553791887,
    0.02017333553791887, 0.02017333553791887, 0.02017333553791887, 0.02017333553791887,
    0.02017333553791887, 0.02017333553791887, 0.02017333553791887, 0.02017333553791887,
    0.02017333553791887, 0.02017333553791887, 0.02017333553791887, 0.02017333553791887,
    0.02017333553791887, 0.02017333553791887, 0.02017333553791887, 0.02017333553791887,
    0.02017333553791887, 0.02017333553791887, 0.02017333553791887])


# For every element a few number are needed to perform an approximate Hirshfeld
# partitioning and to define radial integration grids:
#
# Attributes
# ----------
# populations : np.ndarray, shape=(nshell,)
#
class ElementSetup(object):
    """Parameters that define all molecular grid parameters for one element."""
    def __init__(self, populations, alphas, a, b, npoint):
        r"""Initialize an ElementSetup.

        The electron density of every element is approximately represented by a
        superposition of s-type Slater density functions. One function per shell is added,
        eaching have a population and an exponent. The approximate atomic density has
        spherical symmtery with the following radial dependencies:

        .. math::

            rho(r) = \\sum_{i=0}^{N_\text{shell}} \frac{N_i \alpha_i^3}{8\pi}
                                                   \exp(-\alpha_i r)

        The population and exponent parameters in this module are fitted to densities from
        atomic relativistic LDA calculations, using the MBIS method.

        The radial grid of every element has the following form:

        .. math::

            r(t) = a*t/(1 - b*t)

        where t goes from 0 to npoint-1. a and b are parameters tuned for every element.
        The numbers in this module are borrowed from GPAW.

        Parameters
        ----------
        populations : np.ndarray, shape=(nshell,)
            The populations of the shells, :math:`N_i`
        alphas : np.ndarray, shape=(nshell,)
            The exponents of the shells, :math:`\alpha_i`
        a : float
            The a parameter defined above.
        b : float
            The b parameter defined above.
        npoint : int
            The number of radial integration grid points.
        """
        self.populations = np.asarray(populations)
        self.alphas = np.asarray(alphas)
        self.a = a
        self.b = b
        self.npoint = npoint

    def get_radius_weight(self, t):
        """Return the radius and jacobian for index t."""
        r = self.a*t/(1 - self.b*t)
        w = 4*np.pi*r**2*self.a/(1 - self.b*t)**2
        return r, w

    def compute_density(self, distances):
        result = 0.0
        for ishell in xrange(len(self.populations)):
            prefactor = self.populations[ishell]*self.alphas[ishell]**3/(8*np.pi)
            result += prefactor*np.exp(-distances*self.alphas[ishell])
        return result


element_setups = {
    1: ElementSetup([1.0], [2.0], 0.4/150.0, 1.0/150.0, 150),
    6: ElementSetup([1.70730, 4.29270], [12.79758, 1.85580], 0.4/300.0, 1.0/300.0, 300),
    7: ElementSetup([1.68283, 5.31717], [15.13096, 2.19942], 0.4/300.0, 1.0/300.0, 300),
    8: ElementSetup([1.66122, 6.33878], [17.46129, 2.54326], 0.4/300.0, 1.0/300.0, 300),
}


def generate_molecular_grid(numbers, coordinates):
    """Generate a molecular integration grid.

    Parameters
    ----------
    numbers : np.ndarray, shape=(natom,), dtype=int
        Atomic numbers.
    coordinates : np.ndarray, shape(natom, 3), dtype=float
        Atomic cooorindates.

    Returns
    -------
    points : np.ndarray, shape=(npoint, 3), dtype=float
        Positions of the grid points.
    weights : np.ndarray, shape=(npoint,), dypte=float
        Integration grid weights for all points.
    """
    natom = len(numbers)
    atomic_grids = []
    for iatom0 in xrange(natom):
        # 1) Generate atomic grid points and weights
        setup0 = element_setups[numbers[iatom0]]
        nsphere = len(lebedev_50_weights)
        atomic_points = np.zeros((setup0.npoint, nsphere, 3))
        atomic_weights = np.zeros((setup0.npoint, nsphere))
        for t in xrange(setup0.npoint):
            r, w = setup0.get_radius_weight(t)
            atomic_points[t] = lebedev_50_points*r + coordinates[iatom0]
            atomic_weights[t] = lebedev_50_weights*w
        atomic_points.shape = (-1, 3)
        atomic_weights.shape = (-1,)

        # 2) Evaluate pro-molecular density on the grid. Also keep the density of iatom0.
        pro_mol = 0.0
        pro_atom0 = None
        for iatom1 in xrange(natom):
            deltas = atomic_points - coordinates[iatom1]
            distances = np.sqrt(np.einsum('ij,ij->i', deltas, deltas))
            setup1 = element_setups[numbers[iatom1]]
            pro_atom1 = setup1.compute_density(distances)
            if iatom0 == iatom1:
                pro_atom0 = pro_atom1
            pro_mol += pro_atom1

        # 3) Multiply in the Hirshfeld weight. Simple trick to avoid division by zero.
        ratios = (pro_atom0 + 1e-100)/(pro_mol + 1e-100)
        atomic_weights *= ratios

        # 4) Store
        atomic_grids.append((atomic_points, atomic_weights))

    points = np.concatenate([atomic_points for atomic_points, atomic_weights in atomic_grids])
    weights = np.concatenate([atomic_weights for atomic_points, atomic_weights in atomic_grids])
    return points, weights


def integrate(*args):
    """Replaces grid.integrate for gbasis tests.
       Simply takes a dot product of all 1D numpy arrays passed to it.
    """
    return np.sum(reduce(np.multiply, args))