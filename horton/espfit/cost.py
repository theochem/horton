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
'''ESP cost functions for estimating and testing charges'''


import numpy as np

from horton.log import log
from horton.units import angstrom
from horton.grid.cext import UniformGrid
from horton.espfit.cext import setup_esp_cost_cube, multiply_dens_mask, \
    multiply_near_mask, multiply_far_mask
from horton.utils import typecheck_geo


__all__ = ['ESPCost', 'setup_weights']


class ESPCost(object):
    def __init__(self, A, B, C, natom):
        # Set attributes
        self._A = A
        self._B = B
        self._C = C
        self.natom = natom
        # Rescale parameters not related to atomic charges

    @classmethod
    def from_hdf5(cls, grp):
        return cls(
            grp['A'][:],
            grp['B'][:],
            grp['C'][()],
            grp['natom'][()],
        )

    def to_hdf5(self, grp):
        grp.attrs['class'] = self.__class__.__name__
        grp['A'] = self._A
        grp['B'] = self._B
        grp['C'] = self._C
        grp['natom'] = self.natom

    @classmethod
    def from_grid_data(cls, coordinates, ugrid, vref, weights, rcut=20.0, alpha=None, gcut=None):
        if len(coordinates.shape) != 2 or coordinates.shape[1] != 3:
            raise TypeError('The argument coordinates must be an array with three columns.')
        natom = coordinates.shape[0]
        if alpha is None:
            alpha = 3.0 / rcut
        if gcut is None:
            gcut = 1.1 * alpha
        if isinstance(ugrid, UniformGrid):
            natom = len(coordinates)
            if (ugrid.pbc == [1, 1, 1]).all():
                A = np.zeros((natom+1, natom+1), float)
                B = np.zeros(natom+1, float)
                C = np.zeros((), float)
                setup_esp_cost_cube(ugrid, vref, weights, coordinates, A, B, C, rcut, alpha, gcut)
                return cls(A, B, C, natom)
            else:
                A = np.zeros((natom, natom), float)
                B = np.zeros(natom, float)
                C = np.zeros((), float)
                setup_esp_cost_cube(ugrid, vref, weights, coordinates, A, B, C, 0.0, 0.0, 0.0)
                return cls(A, B, C, natom)
        else:
            raise NotImplementedError

    def value(self, x):
        return np.dot(x, np.dot(self._A, x) - 2*self._B) + self._C

    def value_charges(self, charges):
        if self.natom < len(self._A):
            # Set up a system of equations where all charges are fixed and the
            # remaining parameters are solved for.
            A = self._A[self.natom:,self.natom:]
            B = self._B[self.natom:] - np.dot(charges, self._A[:self.natom,self.natom:])
            C = self._C \
                + np.dot(np.dot(charges, self._A[:self.natom,:self.natom]), charges) \
                - 2*np.dot(self._B[:self.natom], charges)
            x = np.linalg.solve(A, B)
            return C - np.dot(B, x)
        else:
            return self.value(charges)

    def gradient(self, x):
        return 2*(np.dot(self._A, x) - self._B)

    def worst(self, qtot=0.0):
        '''Return a worst-case value for the cost function

           **Optional arguments:**

           qtot
                The total charge of the molecule/crystal

           Higher values for the cost function are still possible but if that
           happens, it is better not to use charges at all.
        '''
        charges = np.zeros(self.natom)
        charges[:] = qtot/self.natom
        return self.value_charges(charges)

    def solve(self, qtot=None, ridge=0.0):
        # apply regularization to atomic degrees of freedom
        A = self._A.copy()
        A.ravel()[::len(A)+1][:self.natom] += ridge*np.diag(A)[:self.natom].mean()
        # construct preconditioned matrices
        norms = np.diag(A)**0.5
        A = A/norms/norms.reshape(-1,1)
        B = self._B/norms

        x = np.linalg.solve(A, B)
        if qtot is not None:
            # Fix the total charge with a lagrange multiplier
            d = np.zeros(len(A))
            d[:self.natom] = 1/norms[:self.natom]
            d[self.natom:] = 0.0
            aid = np.linalg.solve(A, d)
            lagrange = (np.dot(aid, B) - qtot)/np.dot(aid, d)
            x -= aid*lagrange
        x /= norms
        return x


def setup_weights(coordinates, numbers, grid, dens=None, near=None, far=None):
    '''Define a weight function for the ESPCost

       **Arguments:**

       coordinates
            An array with shape (N, 3) containing atomic coordinates.

       numbers
            A vector with shape (N,) containing atomic numbers.

       grid
            A UniformGrid object.

       **Optional arguments:**

       dens
            The density-based criterion. This is a three-tuple with rho, lnrho0
            and sigma. rho is the atomic or the pro-atomic electron density on
            the same grid as the ESP data. lnrho0 and sigma are parameters
            defined in JCTC, 3, 1004 (2007), DOI:10.1021/ct600295n. The weight
            function takes the form::

                exp(-sigma*(ln(rho) - lnrho0)**2)

            Note that the density, rho, should not contain depletions in the
            atomic cores, as is often encountered with pseudo-potential
            computations. In that case it is recommended to construct a
            promolecular density as input for this option.

       near
            Exclude points near the nuclei. This is a dictionary with as items
            (number, (R0, gamma)).

       far
            Exclude points far away. This is a two-tuple: (R0, gamma).
    '''
    natom, coordinates, numbers = typecheck_geo(coordinates, numbers, need_pseudo_numbers=False)
    weights = np.ones(grid.shape)

    # combine three possible mask functions
    if dens is not None:
        log.cite('hu2007', 'for the ESP fitting weight function')
        rho, lnrho0, sigma = dens
        assert (rho.shape == grid.shape).all()
        multiply_dens_mask(rho, lnrho0, sigma, weights)
    if near is not None:
        for i in xrange(natom):
            pair = near.get(numbers[i])
            if pair is None:
                pair = near.get(0)
            if pair is None:
                continue
            r0, gamma = pair
            if r0 > 5*angstrom:
                raise ValueError('The wnear radius is excessive. Please keep it below 5 angstrom.')
            multiply_near_mask(coordinates[i], grid, r0, gamma, weights)
    if far is not None:
        r0, gamma = far
        multiply_far_mask(coordinates, grid, r0, gamma, weights)

    # double that weight goes to zero at non-periodic edges
    return weights
