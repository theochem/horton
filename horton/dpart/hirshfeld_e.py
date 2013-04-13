# -*- coding: utf-8 -*-
# Horton is a Density Functional Theory program.
# Copyright (C) 2011-2012 Toon Verstraelen <Toon.Verstraelen@UGent.be>
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

from horton.cache import just_once
from horton.dpart.base import DPart
from horton.dpart.hirshfeld_i import HirshfeldIDPart
from horton.dpart.linalg import solve_positive
from horton.grid.int1d import SimpsonIntegrator1D
from horton.grid.cext import CubicSpline, dot_multi
from horton.log import log


__all__ = ['HirshfeldEDPart']


# TODO: proofread and add tests for pseudo densities


class HirshfeldEDPart(HirshfeldIDPart):
    '''Extended Hirshfeld partitioning'''
    # TODO: use HEBasis

    name = 'he'
    options = ['local', 'threshold', 'maxiter']

    def _get_proatom_fn(self, index, number, target_charge, first, grid):
        if first:
            return self._proatomdb.get_spline(number)
        else:
            rtf = self._proatomdb.get_rtransform(number)
            int1d = SimpsonIntegrator1D()
            radii = rtf.get_radii()
            weights = (4*np.pi) * radii**2 * int1d.get_weights(len(radii)) * rtf.get_volume_elements()
            assert (weights > 0).all()

            if self.cache.has('records', number):
                records = self.cache.load('records', number)
                pops = self.cache.load('pops', number)
                neq = len(pops)
                A = self.cache.load('A', number)
            else:
                # Construct the basis functions
                basis = []
                #norms = []
                for j0, charge in enumerate(self._proatomdb.get_charges(number, True)):
                    if j0 == 0:
                        # TODO: the corresponding coefficient should be constrained unless pseudo_population==1
                        record = self._proatomdb.get_record(number, charge)
                        radrho = record.rho/record.pseudo_population
                    else:
                        radrho = self._proatomdb.get_record(number, charge).rho - \
                                self._proatomdb.get_record(number, charge+1).rho

                    redundant = False
                    for j1 in xrange(len(basis)):
                        delta = radrho - basis[j1]
                        if dot_multi(weights, delta, delta) < 1e-5:
                            redundant = True
                            break

                    if not redundant:
                        basis.append(radrho)
                        #if j0 == 0:
                        #    norms.append(record.pseudo_population)
                        #else:
                        #    norms.append(1.0)
                basis = np.array(basis)
                #norms = np.array(norms)
                self.cache.dump('basis', number, basis)
                #self.cache.dump('norms', number, norms)

                # Set up system of linear equations:
                neq = len(basis)
                A = np.zeros((neq, neq), float)
                for j0 in xrange(neq):
                    for j1 in xrange(j0+1):
                        A[j0, j1] = dot_multi(weights, basis[j0], basis[j1])
                        A[j1, j0] = A[j0, j1]

                if (np.diag(A) < 0).any():
                    raise ValueError('The diagonal of A must be positive.')

                self.cache.dump('A', number, A)

            B = np.zeros(neq, float)
            tmp = np.zeros(grid.size)
            dens = self.cache.load('mol_dens', index)
            at_weights = self.cache.load('at_weights', index)
            for j0 in xrange(neq):
                tmp[:] = 0.0
                spline = CubicSpline(basis[j0], rtf=rtf)
                grid.eval_spline(spline, self._system.coordinates[index], tmp)
                B[j0] = grid.integrate(at_weights, dens, tmp)

            # Precondition
            scales = np.diag(A)**0.5
            Ap = A/scales/scales.reshape(-1, 1)
            Bp = B/scales

            # Find solution
            target_pseudo_population = self.system.pseudo_numbers[index] - target_charge
            lc_pop = (np.ones(neq)/scales, target_pseudo_population)
            coeffs = solve_positive(Ap, Bp, [lc_pop])
            # correct for scales
            coeffs /= scales

            # Screen output
            if log.do_medium:
                log('                   %10i:&%s' % (index, ' '.join('% 6.3f' % c for c in coeffs)))

            # Construct pro-atom
            proradrho = 0
            for j0 in xrange(neq):
                proradrho += coeffs[j0]*basis[j0]
            error = np.dot(proradrho, weights) - target_pseudo_population
            assert abs(error) < 1e-5

            # Check for negative parts
            if proradrho.min() < 0:
                proradrho[proradrho<0] = 0.0
                error = np.dot(proradrho, weights) - target_pseudo_population
                if log.do_medium:
                    log('                    Pro-atom not positive everywhere. Lost %.5f electrons' % error)

            # Done
            return CubicSpline(proradrho, rtf=rtf)
