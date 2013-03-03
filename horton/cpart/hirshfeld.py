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
from horton.grid.cext import CubicSpline
from horton.dpart.linalg import solve_positive, quadratic_solver


__all__ = [
    'HirshfeldCPart', 'HirshfeldICPart', 'HirshfeldECPart',
]


# TODO: reduce duplicate code


class HirshfeldCPart(CPart):
    name = 'h'

    def __init__(self, system, ui_grid, moldens, proatomdb, scratch, smooth=False):
        '''
           See CPart base class for the description of the arguments.
        '''
        self._proatomdb = proatomdb
        CPart.__init__(self, system, ui_grid, moldens, scratch, smooth)

    def _get_proatomdb(self):
        return self._proatomdb

    proatomdb = property(_get_proatomdb)

    @just_once
    def _init_weight_corrections(self):
        funcs = []
        for i in xrange(self._system.natom):
            n = self._system.numbers[i]

            # TODO: move this check to __init__
            pn = self._system.pseudo_numbers[i]
            pn_expected = self._proatomdb.get_record(n, 0).pseudo_number
            if self._proatomdb.get_record(n, 0).pseudo_number != pn:
                raise ValueError('The pseudo number of atom %i does not match with the proatom database (%i!=%i)' % (i, pn, pn_expected))

            funcs.append((
                self._system.coordinates[i],
                [self._proatomdb.get_spline(n)],
            ))
        wcor = self._ui_grid.compute_weight_corrections(funcs)
        self._scratch.dump('wcor', wcor)

    def _compute_proatom(self, index):
        if log.do_medium:
            log('Computing proatom %i (%i)' % (index, self.system.natom))
        proatom = self._zeros()
        # Construct the pro-atom
        number = self._system.numbers[index]
        center = self._system.coordinates[index]
        spline = self._proatomdb.get_spline(number)
        self._ui_grid.eval_spline(spline, center, proatom)
        proatom += 1e-100
        self._scratch.dump('at_weights', index, proatom)

    def _compute_promoldens(self):
        promoldens = self._zeros()
        proatom = self._zeros()
        for i in xrange(self._system.natom):
            self._scratch.load('at_weights', i, output=proatom)
            promoldens += proatom
        self._scratch.dump('promoldens', promoldens)

    def _compute_at_weights(self):
        promoldens = self._scratch.load('promoldens')
        proatom = self._zeros()
        for i in xrange(self._system.natom):
            self._scratch.load('at_weights', i, output=proatom)
            proatom /= promoldens
            self._scratch.dump('at_weights', i, proatom)

    @just_once
    def _init_at_weights(self):
        # A) Compute all the pro-atomic densities.
        for i in xrange(self._system.natom):
            self._compute_proatom(i)
        self._compute_promoldens()
        self._compute_at_weights()

    def do_all(self):
        '''Computes all reasonable properties and returns a corresponding list of keys'''
        names = CPart.do_all(self)
        self.do_dispersion()
        return names + ['volume_ratios', 'c6s']

    @just_once
    def do_dispersion(self):
        # Method by Alexandre Tkatchenko and Matthias Scheffler:
        #   PhysRevLett-v102-p073005-y2009.pdf
        #   http://dx.doi.org/10.1103/PhysRevLett.102.073005
        # Reference C6 values taken from X. Chu and A. Dalgarno:
        #   JChemPhys-v121-p4083-y2004.pdf
        #   http://dx.doi.org/10.1063/1.1779576
        #   (Corrected values from tables VII to XII)
        # C6 value for Hydrogen taken from Zong-Chao Yan, James F. Babb, A. Dalgarno and G. W. F. Drake
        #   PhysRevA-v54-p2824-y1996.pdf
        #   http://dx.doi.org/10.1103/PhysRevA.54.2824
        ref_c6s = { # reference C6 values in atomic units
            1: 6.499, 2: 1.42, 3: 1392.0, 4: 227.0, 5: 99.5, 6: 46.6, 7: 24.2,
            8: 15.6, 9: 9.52, 10: 6.20, 11: 1518.0, 12: 626.0, 13: 528.0, 14:
            305.0, 15: 185.0, 16: 134.0, 17: 94.6, 18: 64.2, 19: 3923.0, 20:
            2163.0, 21: 1383.0, 22: 1044.0, 23: 832.0, 24: 602.0, 25: 552.0, 26:
            482.0, 27: 408.0, 28: 373.0, 29: 253.0, 30: 284.0, 31: 498.0, 32:
            354.0, 33: 246.0, 34: 210.0, 35: 162.0, 36: 130.0, 37: 4769.0, 38:
            3175.0, 49: 779.0, 50: 659.0, 51: 492.0, 52: 445.0, 53: 385.0,
        }

        volume_ratios, new_volume_ratios = self._cache.load('volume_ratios', alloc=self.system.natom)
        c6s, new_c6s = self._cache.load('c6s', alloc=self.system.natom)

        if new_volume_ratios or new_c6s:
            self.do_volumes()
            volumes = self._cache.load('volumes')

            if log.do_medium:
                log('Computing atomic dispersion coefficients.')

            for i in xrange(self.system.natom):
                n = self.system.numbers[i]
                if n not in ref_c6s:
                    raise NotImplementedError('No reference C6 value availlable for atom number %i.' % n)
                ref_volume = self.proatomdb.get_record(n, 0).get_moment(3)/n
                volume_ratios[i] = volumes[i]/ref_volume
                c6s[i] = (volume_ratios[i])**2*ref_c6s[n]


class HirshfeldICPart(HirshfeldCPart):
    name = 'hi'
    options = ['smooth', 'max_iter', 'threshold']

    def __init__(self, system, ui_grid, moldens, proatomdb, scratch, smooth=False, max_iter=100, threshold=1e-4):
        '''
           **Optional arguments:** (those not present in the base class)

           max_iter
                The maximum number of iterations. If no convergence is reached
                in the end, no warning is given.

           threshold
                The procedure is considered to be converged when the maximum
                change of the charges between two iterations drops below this
                threshold.
        '''
        self._max_iter = max_iter
        self._threshold = threshold
        HirshfeldCPart.__init__(self, system, ui_grid, moldens, proatomdb, scratch, smooth)

    def _get_isolated_atom(self, i, charge, output):
        key = ('isolated_atom', i, charge)
        if key in self._scratch:
            self._scratch.load(*key, output=output)
        else:
            number = self._system.numbers[i]
            center = self._system.coordinates[i]
            spline = self._proatomdb.get_spline(number, charge)
            if log.do_medium:
                log('Computing isolated atom %i (n=%i, q=%+i)' % (i, number, charge))
            output[:] = 0
            self._ui_grid.eval_spline(spline, center, output)
            self._scratch.dump(*key, data=output)

    def _compute_proatom(self, i, target_charge):
        # Pro-atoms are (temporarily) stored in at_weights for efficiency.
        proatom = self._zeros()
        # Construct the pro-atom
        icharge = int(np.floor(target_charge))
        self._get_isolated_atom(i, icharge, proatom)
        if float(icharge) != target_charge:
            x = target_charge - icharge
            proatom *= 1-x
            ipseudo_pop = self.system.pseudo_numbers[i] - icharge
            if ipseudo_pop > 1:
                tmp = self._zeros()
                self._get_isolated_atom(i, icharge+1, tmp)
                tmp *= x
                proatom += tmp
        proatom += 1e-100 # avoid division by zero
        self._scratch.dump('at_weights', i, proatom)

    @just_once
    def _init_at_weights(self):
        old_charges = np.zeros(self._system.natom)
        if log.medium:
            log.hline()
            log('Iteration       Change')
            log.hline()

        counter = 0
        change = 1.0
        while True:
            counter += 1

            # Update pro-atoms and populations
            for i in xrange(self._system.natom):
                self._compute_proatom(i, old_charges[i])
            self._compute_promoldens()
            self._compute_at_weights()
            new_populations = self._compute_populations()
            new_charges = self.system.numbers - new_populations

            # Check for convergence
            change = abs(new_charges - old_charges).max()
            if log.medium:
                log('%9i   %10.5e' % (counter, change))
            if change < self._threshold or counter >= self._max_iter:
                break

            old_charges = new_charges

        if log.medium:
            log.hline()

        self._cache.dump('populations', new_populations)
        self._cache.dump('niter', counter)
        self._cache.dump('change', change)

    def do_all(self):
        names = HirshfeldCPart.do_all(self)
        return names + ['niter', 'change']


class HirshfeldECPart(HirshfeldICPart):
    name = 'he'

    @just_once
    def _init_weight_corrections(self):
        HirshfeldICPart._init_weight_corrections(self)

        funcs = []
        for i in xrange(self._system.natom):
            center = self._system.coordinates[i]
            splines = []
            number = self._system.numbers[i]
            charges = self._proatomdb.get_charges(number, safe=True)
            rtf = self._proatomdb.get_rtransform(number)
            splines = []
            for j0 in xrange(len(charges)-1):
                rad0 = self._proatomdb.get_record(number, charges[j0+1]).rho - \
                       self._proatomdb.get_record(number, charges[j0]).rho
                for j1 in xrange(j0+1):
                    rad1 = self._proatomdb.get_record(number, charges[j1+1]).rho - \
                           self._proatomdb.get_record(number, charges[j1]).rho
                    splines.append(CubicSpline(rad0*rad1, rtf=rtf))
            funcs.append((center, splines))
        wcor_fit = self._ui_grid.compute_weight_corrections(funcs)
        self._scratch.dump('wcor_fit', wcor_fit)

    def _get_constant_fn(self, i, output):
        key = ('isolated_atom', i, 0)
        if key in self._scratch:
            self._scratch.load(*key, output=output)
        else:
            number = self._system.numbers[i]
            center = self._system.coordinates[i]
            spline = self._proatomdb.get_spline(number)
            if log.do_medium:
                log('Computing constant fn %i (n=%i)' % (i, number))
            output[:] = 0
            self._ui_grid.eval_spline(spline, center, output)
            self._scratch.dump(*key, data=output)

    def _get_basis_fn(self, i, j, output):
        # A local cache is used that only exists during _init_at_weights:
        key = ('basis', i, j)
        if key in self._scratch:
            self._scratch.load(*key, output=output)
        else:
            number = self._system.numbers[i]
            charges = self._proatomdb.get_charges(number, safe=True)
            center = self._system.coordinates[i]
            spline = self._proatomdb.get_spline(number, {charges[j+1]: 1, charges[j]: -1})
            if log.do_medium:
                log('Computing basis fn %i_%i (n=%i, q: %+i_%+i)' % (i, j, number, charges[j+1], charges[j]))
            output[:] = 0
            self._ui_grid.eval_spline(spline, center, output)
            self._scratch.dump(*key, data=output)

    def _compute_proatom(self, i, target_charge):
        # Pro-atoms are (temporarily) stored in at_weights for efficiency.
        if ('at_weights', i) not in self._scratch:
            # just use the Hirshfeld definition to get started
            proatom = self._zeros()
            self._get_constant_fn(i, proatom)
            proatom += 1e-100
            self._scratch.dump('at_weights', i, proatom)
        else:
            # 1) construct aim in temporary array
            at_weights = self._scratch.load('at_weights', i)
            aimdens = self._scratch.load('moldens')
            aimdens *= at_weights

            #    rename variables for clarity
            work0 = at_weights
            del at_weights
            work1 = self._zeros()

            # 2) setup equations
            number = self._system.numbers[i]
            charges = np.array(self._proatomdb.get_charges(number, safe=True))
            neq = len(charges)-1

            #    get the weight corrections for the least squares problem
            # TODO: default value for load
            if 'wcor_fit' in self._scratch:
                wcor_fit = self._scratch.load('wcor_fit')
            else:
                wcor_fit = None

            #    compute A
            A, new = self._cache.load('A', i, alloc=(neq, neq))
            if new:
                for j0 in xrange(neq):
                    self._get_basis_fn(i, j0, work0)
                    for j1 in xrange(j0+1):
                        self._get_basis_fn(i, j1, work1)
                        A[j0,j1] = self._ui_grid.integrate(work0, work1, wcor_fit)
                        A[j1,j0] = A[j0,j1]


                #    precondition the equations
                scales = np.diag(A)**0.5
                A /= scales
                A /= scales.reshape(-1, 1)
                if log.do_medium:
                    evals = np.linalg.eigvalsh(A)
                    cn = abs(evals).max()/abs(evals).min()
                    sr = abs(scales).max()/abs(scales).min()
                    log('                   %10i: CN=%.5e SR=%.5e' % (i, cn, sr))
                self._cache.dump('scales', i, scales)
            else:
                scales = self._cache.load('scales', i)

            #    compute B and precondition
            self._get_constant_fn(i, work1)
            aimdens -= work1
            B = np.zeros(neq, float)
            for j0 in xrange(neq):
                self._get_basis_fn(i, j0, work0)
                B[j0] = self._ui_grid.integrate(aimdens, work0, wcor_fit)
            B /= scales

            # 3) find solution
            lc_pop = (np.ones(neq)/scales, -target_charge)
            #coeffs = solve_positive(A, B, [lc_pop])
            lcs_par = []
            for j0 in xrange(neq):
                lc = np.zeros(neq)
                lc[j0] = 1.0/scales[j0]
                lcs_par.append((lc, -1))
            coeffs = quadratic_solver(A, B, [lc_pop], lcs_par, rcond=0)
            # correct for scales
            coeffs /= scales

            if log.do_medium:
                log('                   %10i:&%s' % (i, ' '.join('% 6.3f' % c for c in coeffs)))

            # 4) construct the pro-atom
            proatom = work1
            del work1

            for j0 in xrange(neq):
                work0[:] = 0
                self._get_basis_fn(i, j0, work0)
                work0 *= coeffs[j0]
                proatom += work0

            proatom += 1e-100
            self._scratch.dump('at_weights', i, proatom)
