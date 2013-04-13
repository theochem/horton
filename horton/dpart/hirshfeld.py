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
from horton.dpart.linalg import solve_positive
from horton.grid.int1d import SimpsonIntegrator1D
from horton.grid.cext import CubicSpline, dot_multi
from horton.log import log


__all__ = ['HirshfeldDPart', 'HirshfeldIDPart', 'HirshfeldEDPart']


# TODO: proofread and add tests for pseudo densities


class HirshfeldDPart(DPart):
    '''Base class for Hirshfeld partitioning'''
    def __init__(self, molgrid, proatomdb, local=True):
        self._proatomdb = proatomdb
        self._pro_mol_valid = False
        DPart.__init__(self, molgrid, local)

    def _get_proatomdb(self):
        return self._proatomdb

    proatomdb = property(_get_proatomdb)

    def _init_log(self):
        DPart._init_log(self)
        if log.do_medium:
            log.deflist([
                ('Scheme', 'Hirshfeld'),
                ('Proatomic DB',  self._proatomdb),
            ])
            log.cite('hirshfeld1977', 'the use of Hirshfeld partitioning')

    def _at_weights_helper(self, i0, grid, at_weights):
        # Load work arrays from cache for efficiency
        (work, pro_mol), new = self.cache.load('at_weights_work', grid.size, alloc=(2, grid.size))
        if not new:
            pro_mol[:] = 0.0
        if self.local or not self._pro_mol_valid:
            # In case of local grids, the pro-molecule must be recomputed for
            # every grid. In case of a global grid, this only happens the
            # first time after the pro-atoms were updated.
            for i1 in xrange(self.system.natom):
                proatom_fn = self.cache.load('proatom_fn', i1)
                work[:] = 0.0
                grid.eval_spline(proatom_fn, self.system.coordinates[i1], work)
                if i1 == i0:
                    at_weights[:] = work
                pro_mol[:] += work
            # The following seems worse than it is. It does nothing to the
            # relevant numbers. It just avoids troubles in the division.
            pro_mol[:] += 1e-100
        else:
            # In case of a global grid, and when the pro-molecule is up to date,
            # only the pro-atom needs to be recomputed.
            proatom_fn = self.cache.load('proatom_fn', i0)
            at_weights[:] = 0.0
            grid.eval_spline(proatom_fn, self.system.coordinates[i0], at_weights)
        # Finally compute the ratio
        at_weights[:] /= pro_mol

    def _at_weights_cleanup(self):
        # Get rid of cached work arrays
        for i, grid in self.iter_grids():
            self.cache.discard('at_weights_work', id(grid.size))

    @just_once
    def _init_at_weights(self):
        for i in xrange(self.system.natom):
            proatom_fn = self.cache.load('proatom_fn', i, default=None)
            if proatom_fn is None:
                n = self.system.numbers[i]
                proatom_fn = self._proatomdb.get_spline(n)
                self.cache.dump('proatom_fn', i, proatom_fn)

        for i0, grid in self.iter_grids():
            at_weights, new = self.cache.load('at_weights', i0, alloc=grid.size)
            if new:
                self._at_weights_helper(i0, grid, at_weights)

        self._at_weights_cleanup()


class HirshfeldIDPart(HirshfeldDPart):
    '''Iterative Hirshfeld partitioning'''
    def __init__(self, molgrid, proatomdb, local=True, threshold=1e-4, maxiter=500):
        self._threshold = threshold
        self._maxiter = maxiter
        HirshfeldDPart.__init__(self, molgrid, proatomdb, local)

    def _init_log(self):
        DPart._init_log(self)
        if log.do_medium:
            log.deflist([
                ('Scheme', 'Hirshfeld-I'),
                ('Convergence threshold', '%.1e' % self._threshold),
                ('Maximum iterations', self._maxiter),
                ('Proatomic DB',  self._proatomdb),
            ])
            log.cite('bultinck2007', 'the use of Hirshfeld-I partitioning')

    def _proatom_change(self, number, old, new):
        '''Compute the difference between an old and a new spline

           **Arguments:**

           number
                Element number

           old, new
                The old and new splines
        '''
        rtf = self.proatomdb.get_rtransform(number)
        return (4*np.pi)*dot_multi(
            rtf.get_radii()**2, # TODO: get routines are slow
            rtf.get_volume_elements(),
            SimpsonIntegrator1D().get_weights(rtf.npoint),
            (old.copy_y() - new.copy_y())**2 # TODO: copy is slow
        )

    def _get_proatom_fn(self, index, number, target_charge, first, grid):
        icharge = int(np.floor(target_charge))
        x = target_charge - icharge
        # check if icharge record should exist
        pseudo_pop = self.system.pseudo_numbers[index] - icharge
        if pseudo_pop == 1:
            return self._proatomdb.get_spline(number, {icharge: 1-x})
        elif pseudo_pop > 1:
            return self._proatomdb.get_spline(number, {icharge: 1-x, icharge+1: x})
        else:
            raise ValueError('Requesting a pro-atom with a negative (pseudo) population')

    @just_once
    def _init_at_weights(self):
        # Perform one general check in the beginning to keep things simple.
        if all(self.cache.has('at_weights', i) for i in xrange(self.system.natom)):
            return

        self.do_mol_dens()

        # Iterative Hirshfeld loop
        charges = np.zeros(self.system.natom)
        pseudo_populations = np.zeros(self.system.natom)
        counter = 0
        if log.do_medium:
            log('Iterative Hirshfeld partitioning loop')
            log.hline()
            log('Counter      Change   QTot Error')
            log.hline()

        while True:
            # Update current pro-atoms
            change = 0.0
            first_iter = True
            for i, grid in self.iter_grids():
                old_proatom_fn = self.cache.load('proatom_fn', i, default=None)
                number = self.system.numbers[i]
                charge = charges[i]
                proatom_fn = self._get_proatom_fn(i, number, charge, old_proatom_fn is None, grid)

                self.cache.dump('proatom_fn', i, proatom_fn)
                if old_proatom_fn is not None:
                    first_iter = False
                    change += self._proatom_change(number, old_proatom_fn, proatom_fn)

            # Enforce (single) update of pro-molecule in case of a global grid
            if not self.local:
                self._pro_mol_valid = False

            # Compute populations
            for i, grid in self.iter_grids():
                # Compute weight
                at_weights, new = self.cache.load('at_weights', i, alloc=grid.size)
                self._at_weights_helper(i, grid, at_weights)

                # Compute population
                dens = self.cache.load('mol_dens', i)
                pseudo_populations[i] = grid.integrate(at_weights, dens)
            charges = self.system.pseudo_numbers - pseudo_populations

            change = np.sqrt(change/self.system.natom)
            if log.do_medium:
                pop_error = charges.sum() - self.system.charge
                log('%7i  %10.5e  %10.5e' % (counter, change, pop_error))

            if not first_iter and change < self._threshold:
                self._converged = True
                break

            if counter > self._maxiter:
                break

            counter += 1

        if log.do_medium:
            log.hline()
            log('Converged: %s' % self._converged)
            log.blank()

        self._at_weights_cleanup()
        populations = pseudo_populations - (self.system.pseudo_numbers - self.system.numbers)
        self.cache.dump('pseudo_populations', pseudo_populations)
        self.cache.dump('populations', populations)
        self.cache.dump('charges', charges)


class HirshfeldEDPart(HirshfeldIDPart):
    '''Extended Hirshfeld partitioning'''

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
