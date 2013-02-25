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

from horton.cpart.base import CPart
from horton.cache import just_once, Cache
from horton.log import log
from horton.dpart.linalg import solve_positive, quadratic_solver


__all__ = [
    'HirshfeldCPart', 'HirshfeldICPart', 'HirshfeldECPart',
]


# TODO: maxiter, tolerance parameter
# TODO: reduce duplicate code


class HirshfeldCPart(CPart):
    name = 'h'

    def __init__(self, system, ui_grid, moldens, proatomdb, smooth):
        '''
           **Arguments:**

           system
                The system to be partitioned.

           ui_grid
                The uniform integration grid based on the cube file.

           moldens
                The all-electron molecular density grid data.

           proatomdb
                The database of proatomic densities.

           smooth
                When set to True, no corrections are included to integrate
                the cusps.
        '''
        self._proatomdb = proatomdb
        CPart.__init__(self, system, ui_grid, moldens, smooth)

    def _get_proatomdb(self):
        return self._proatomdb

    proatomdb = property(_get_proatomdb)

    @just_once
    def _init_weight_corrections(self):
        if self.smooth:
            self._cache.dump('wcor', None)
            return

        funcs = []
        for i in xrange(self._system.natom):
            n = self._system.numbers[i]
            pn = self._system.pseudo_numbers[i]
            pn_expected = self._proatomdb.get_record(n, 0).pseudo_number
            if self._proatomdb.get_record(n, 0).pseudo_number != pn:
                raise ValueError('The pseudo number of atom %i does not match with the proatom database (%i!=%i)' % (i, pn, pn_expected))
            funcs.append((
                self._system.coordinates[i],
                [(('isolated_atom', i, 0), self._proatomdb.get_spline(n), float(pn))],
            ))
        self._ui_grid.compute_weight_corrections(funcs, cache=self._cache)

    def _compute_proatom(self, index):
        if self._cache.has('isolated_atom', index, 0):
            self._cache.rename(('isolated_atom', index, 0), ('at_weights', index))
            proatom = self._cache.load('at_weights', index)
        else:
            if log.do_medium:
                log('Computing proatom %i (%i)' % (index, self.system.natom))
            proatom, new = self._cache.load('at_weights', index, alloc=self._ui_grid.shape)
            # Construct the pro-atom
            assert new
            number = self._system.numbers[index]
            center = self._system.coordinates[index]
            spline = self._proatomdb.get_spline(number)
            self._ui_grid.eval_spline(spline, center, proatom)
        proatom += 1e-100

    def _compute_promoldens(self):
        promoldens, new = self._cache.load('promoldens', alloc=self._ui_grid.shape)
        if not new:
            promoldens[:] = 0
        for i in xrange(self._system.natom):
            proatom = self._cache.load('at_weights', i)
            promoldens += proatom

    def _compute_at_weights(self):
        promoldens = self._cache.load('promoldens')
        for i in xrange(self._system.natom):
            proatom = self._cache.load('at_weights', i)
            proatom /= promoldens

    @just_once
    def _init_at_weights(self):
        # A) Compute all the pro-atomic densities.
        for i in xrange(self._system.natom):
            self._compute_proatom(i)
        self._compute_promoldens()
        self._compute_at_weights()
        self._cache.discard('promoldens')



class HirshfeldICPart(HirshfeldCPart):
    name = 'hi'

    def _get_isolated_atom(self, i, charge):
        # A local cache is used that only exists during _init_at_weights:
        isolated_atom, new = self._cache.load('isolated_atom', i, charge, alloc=self._ui_grid.shape)
        if new:
            number = self._system.numbers[i]
            center = self._system.coordinates[i]
            spline = self._proatomdb.get_spline(number, charge)
            if log.do_medium:
                log('Computing isolated atom %i (n=%i, q=%+i)' % (i, number, charge))
            self._ui_grid.eval_spline(spline, center, isolated_atom)
        return isolated_atom

    def _compute_proatom(self, i, target_charge):
        # Pro-atoms are (temporarily) stored in at_weights for efficiency.
        proatom, new = self._cache.load('at_weights', i, alloc=self._ui_grid.shape)
        # Construct the pro-atom
        icharge = int(np.floor(target_charge))
        proatom[:] = self._get_isolated_atom(i, icharge)
        if float(icharge) != target_charge:
            x = target_charge - icharge
            proatom *= 1-x
            ipseudo_pop = self.system.pseudo_numbers[i] - icharge
            if ipseudo_pop > 1:
                tmp = self._local_cache.load('tmp', alloc=self._ui_grid.shape)[0]
                tmp[:] = self._get_isolated_atom(i, icharge+1)
                tmp *= x
                proatom += tmp
        proatom += 1e-100 # avoid division by zero
        return proatom

    @just_once
    def _init_at_weights(self):
        self._local_cache = Cache()

        old_charges = np.zeros(self._system.natom)
        if log.medium:
            log.hline()
            log('Iteration       Change')
            log.hline()
        for counter in xrange(500):
            # Construct pro-atoms
            for i in xrange(self._system.natom):
                self._compute_proatom(i, old_charges[i])
            self._compute_promoldens()
            self._compute_at_weights()
            new_populations = self._compute_populations()
            new_charges = self.system.numbers - new_populations
            change = abs(new_charges - old_charges).max()
            if log.medium:
                log('%9i   %10.5e' % (counter, change))
            if change < 1e-4:
                break

            old_charges = new_charges

        if log.medium:
            log.hline()

        self._cache.dump('populations', new_populations)
        del self._local_cache


class HirshfeldECPart(HirshfeldICPart):
    name = 'he'

    @just_once
    def _init_weight_corrections(self):
        if self.smooth:
            self._cache.dump('wcor', None)
            self._cache.dump('wcor_fit', None)
            return

        HirshfeldICPart._init_weight_corrections(self)

        from horton.grid import SimpsonIntegrator1D, CubicSpline
        int1d = SimpsonIntegrator1D()
        funcs = []
        for i in xrange(self.system.natom):
            n = self.system.numbers[i]
            rtf = self.proatomdb.get_rtransform(n)
            radii = rtf.get_radii()
            weights1d = (4*np.pi) * radii**2 * rtf.get_volume_elements() * int1d.get_weights(len(radii))

            charges = self._proatomdb.get_charges(n)
            rho = self._proatomdb.get_record(n, charges[0]).rho
            spline = CubicSpline(rho**2, rtf=rtf)
            int_exact = np.dot(rho**2, weights1d)

            funcs.append((
                self._system.coordinates[i],
                [(None, spline, int_exact)],
            ))

        wcor_fit = self._ui_grid.compute_weight_corrections(funcs)
        self._cache.dump('wcor_fit', wcor_fit)

    def _get_basis_fn(self, i, j):
        # A local cache is used that only exists during _init_at_weights:
        isolated_atom, new = self._cache.load('basis', i, j, alloc=self._ui_grid.shape)
        number = self._system.numbers[i]
        charges = self._proatomdb.get_charges(number)
        if new:
            number = self._system.numbers[i]
            center = self._system.coordinates[i]
            if j == 0:
                pop = self.proatomdb.get_record(number, charges[0]).pseudo_population
                spline = self._proatomdb.get_spline(number, {charges[0]: 1.0/pop})
            else:
                spline = self._proatomdb.get_spline(number, {charges[j]: 1, charges[j-1]: -1})
            if log.do_medium:
                log('Computing basisfn %i_%i (n=%i, q=%+i)' % (i, j, number, charges[j]))
            self._ui_grid.eval_spline(spline, center, isolated_atom)
        return isolated_atom

    def _compute_proatom(self, i, target_charge):
        # Pro-atoms are (temporarily) stored in at_weights for efficiency.
        if self._cache.has('isolated_atom', i, 0):
            # just use the Hirshfeld definition to get started
            self._cache.rename(('isolated_atom', i, 0), ('at_weights', i))
            proatom = self._cache.load('at_weights', i)
            proatom += 1e-100
        else:
            wcor_fit = self._cache.load('wcor_fit')
            proatom, new = self._cache.load('at_weights', i, alloc=self._ui_grid.shape)
            # do a least-squares fit to the previos AIM

            # 1) construct aim in temporary array
            aimdens = self._local_cache.load('aimdens', alloc=self._ui_grid.shape)[0]
            aimdens[:] = proatom
            aimdens *= self._cache.load('moldens')

            # 2) setup equations
            number = self._system.numbers[i]
            charges = np.array(self._proatomdb.get_charges(number))
            neq = len(charges)

            A, new = self._local_cache.load('fit_A', i, alloc=(neq, neq))
            if new:
                for j0 in xrange(neq):
                    basis0 = self._get_basis_fn(i, j0)
                    for j1 in xrange(j0+1):
                        basis1 = self._get_basis_fn(i, j1)
                        A[j0,j1] = self._ui_grid.integrate(basis0, basis1, wcor_fit)
                        A[j1,j0] = A[j0,j1]
                # precondition the equations
                scales = np.diag(A)**0.5
                A /= scales
                A /= scales.reshape(-1, 1)
                self._local_cache.dump('fit_scales', i, scales)
                if log.do_medium:
                    evals = np.linalg.eigvalsh(A)
                    cn = abs(evals).max()/abs(evals).min()
                    sr = abs(scales).max()/abs(scales).min()
                    log('                   %10i: CN=%.5e SR=%.5e' % (i, cn, sr))
            else:
                scales = self._local_cache.load('fit_scales', i)

            B = np.zeros(neq, float)
            for j0 in xrange(neq):
                basis0 = self._get_basis_fn(i, j0)
                B[j0] = self._ui_grid.integrate(aimdens, basis0, wcor_fit)
            B /= scales


            # 3) find solution
            target_pseudo_pop = self.system.pseudo_numbers[i] - target_charge
            lc_pop = (np.ones(neq)/scales, target_pseudo_pop)
            coeffs = solve_positive(A, B, [lc_pop])
            # correct for scales
            coeffs /= scales

            if log.do_medium:
                log('                   %10i:&%s' % (i, ' '.join('% 6.3f' % c for c in coeffs)))

            # 4) construct the pro-atom
            tmp = aimdens
            del aimdens
            proatom[:] = 0
            for j0 in xrange(neq):
                tmp[:] = self._get_basis_fn(i, j0)
                tmp *= coeffs[j0]
                proatom += tmp
        proatom += 1e-100
