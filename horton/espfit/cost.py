# -*- coding: utf-8 -*-
# Horton is a Density Functional Theory program.
# Copyright (C) 2011-2013 Toon Verstraelen <Toon.Verstraelen@UGent.be>
#
# This file is part of Horton.
#
# Horton is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
#
# Horton is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, see <http://www.gnu.org/licenses/>
#
#--


import numpy as np

from horton.grid.cext import UniformIntGrid
from horton.espfit.cext import setup_esp_cost_cube, multiply_dens_mask, \
    multiply_near_mask, multiply_far_mask


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
    def from_grid_data(cls, system, grid, vref, weights, rcut=20, alpha_scale=3.0, gcut_scale=1.1):
        alpha = alpha_scale / rcut
        gcut = gcut_scale * alpha
        if isinstance(grid, UniformIntGrid):
            A = np.zeros((system.natom+1, system.natom+1), float)
            B = np.zeros(system.natom+1, float)
            C = np.zeros((), float)
            setup_esp_cost_cube(grid, vref, weights, system.coordinates, A, B, C, rcut, alpha, gcut)
            return cls(A, B, C, system.natom)
        else:
            raise NotImplementedError

    def value(self, x):
        return np.dot(x, np.dot(self._A, x) - 2*self._B) + self._C

    def value_charges(self, charges):
        if self.natom < len(self._A):
            # Set up a system where all charges are fixed and the remaning
            # parameters are solved for.
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
                The total charge of the system

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


def setup_weights(system, grid, dens=None, near=None, far=None):
    '''Define a weight function for the ESPCost

       **Arguments:**

       system
            The system for which the weight function must be defined

       grid
            A UniformGrid object.

       **Optional arguments:**

       dens
            The density-based criterion. This is a three-tuple with rho, rho0
            and alpha.

       near
            Exclude points near the nuclei. This is a dictionary with as items
            (number, (R0, gamma)).

       far
            Exclude points far away. This is a two-tuple: (R0, gamma).
    '''
    weights = np.ones(grid.shape)

    # combine three possible mask functions
    if dens is not None:
        rho, rho0, alpha = dens
        assert (rho.shape == grid.shape).all()
        multiply_dens_mask(rho, rho0, alpha, weights)
    if near is not None:
        for i in xrange(system.natom):
            pair = near.get(system.numbers[i])
            if pair is None:
                continue
            r0, gamma = pair
            multiply_near_mask(system.coordinates[i], grid, r0, gamma, weights)
    if far is not None:
        r0, gamma = far
        multiply_far_mask(system.coordinates, grid, r0, gamma, weights)

    # double that weight goes to zero at non-periodic edges
    return weights
