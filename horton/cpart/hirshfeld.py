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

class StockholderCPart(CPart):
    name = 'h'

    def __init__(self, system, ui_grid, moldens, store, smooth=False):
        '''
           See CPart base class for the description of the arguments.
        '''
        CPart.__init__(self, system, ui_grid, moldens, store, smooth)
        assert 'promoldens' in self._cache

    def compute_proatom(self, index, output):
        # This routine must compute the pro-atom, not from scratch, but based
        # on the info available after the partitioning. The default behavior
        # is to load the pro-atom from the store, but in case the store is fake
        # subclasses must be able to reconstruct the pro-atom on the fly.
        raise NotImplementedError

    def compute_at_weights(self, index, output):
        self.compute_proatom(index, output)
        output /= self._cache.load('promoldens')


class HirshfeldCPart(StockholderCPart):
    name = 'h'

    def __init__(self, system, ui_grid, moldens, proatomdb, store, smooth=False):
        '''
           See CPart base class for the description of the arguments.
        '''
        self._proatomdb = proatomdb
        CPart.__init__(self, system, ui_grid, moldens, store, smooth)

    def _get_proatomdb(self):
        return self._proatomdb

    proatomdb = property(_get_proatomdb)

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
        self._cache.dump('wcor', wcor)

    def compute_proatom(self, index, output):
        if log.do_medium:
            log('Computing proatom %i (%i)' % (index, self.system.natom))
        # Construct the pro-atom
        number = self._system.numbers[index]
        center = self._system.coordinates[index]
        spline = self._proatomdb.get_spline(number)
        output[:] = 0.0
        self._ui_grid.eval_spline(spline, center, output)
        output += 1e-100

    def _init_partitioning(self):
        work = self.zeros()
        self._update_promolecule(work)
        self._store_at_weights(work)

    def _update_promolecule(self, work):
        promoldens = self._cache.load('promoldens', alloc=self._ui_grid.shape)[0]
        promoldens[:] = 0.0
        for index in xrange(self._system.natom):
            self.compute_proatom(index, work)
            promoldens += work
            # Merely for efficiency:
            self._store.dump(work, 'at_weights', index)

    def _store_at_weights(self, work):
        # Compute the atomic weight functions if this is useful. This is merely
        # a matter of efficiency.
        promoldens = self._cache.load('promoldens')
        for index in xrange(self._system.natom):
            if ('at_weights', index) in self._store:
                self._store.load(work, 'at_weights', index)
                work /= promoldens
                self._store.dump(work, 'at_weights', index)


class HirshfeldICPart(HirshfeldCPart):
    name = 'hi'
    options = ['smooth', 'max_iter', 'threshold']

    def __init__(self, system, ui_grid, moldens, proatomdb, store, smooth=False, max_iter=100, threshold=1e-4):
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
        HirshfeldCPart.__init__(self, system, ui_grid, moldens, proatomdb, store, smooth)

    def _get_isolated_atom(self, i, charge, output):
        key = ('isolated_atom', i, charge)
        if key in self._store:
            self._store.load(output, *key)
        else:
            number = self._system.numbers[i]
            center = self._system.coordinates[i]
            spline = self._proatomdb.get_spline(number, charge)
            if log.do_medium:
                log('Computing isolated atom %i (n=%i, q=%+i)' % (i, number, charge))
            output[:] = 0.0
            self._ui_grid.eval_spline(spline, center, output)
            self._store.dump(output, *key)

    def compute_proatom(self, i, output):
        # Construct the pro-atom
        target_charge = self._cache.load('charges')[i]
        icharge = int(np.floor(target_charge))
        self._get_isolated_atom(i, icharge, output)
        if float(icharge) != target_charge:
            x = target_charge - icharge
            output *= 1-x
            ipseudo_pop = self.system.pseudo_numbers[i] - icharge
            if ipseudo_pop > 1:
                tmp = self.zeros()
                self._get_isolated_atom(i, icharge+1, tmp)
                tmp *= x
                output += tmp
        output += 1e-100 # avoid division by zero

    def _update_proatom_pars(self):
        pass

    def _init_partitioning(self):
        charges = self._cache.load('charges', alloc=self._system.natom)[0]
        if log.medium:
            log.hline()
            log('Iteration       Change')
            log.hline()

        counter = 0
        change = 1.0
        work = self.zeros()

        while True:
            counter += 1

            # Update the pro-molecule density
            self._update_promolecule(work)

            # Compute the atomic weight functions if this is useful. This is merely
            # a matter of efficiency.
            self._store_at_weights(work)

            # Update the parameters that determine the pro-atoms.
            self._update_proatom_pars()

            # Compute the charges
            old_charges = charges.copy()
            for index in xrange(self._system.natom):
                pseudo_population = self.compute_pseudo_population(index, work)
                charges[index] = self.system.pseudo_numbers[index] - pseudo_population

            # Check for convergence
            change = abs(charges - old_charges).max()
            if log.medium:
                log('%9i   %10.5e' % (counter, change))
            if change < self._threshold or counter >= self._max_iter:
                break

        if log.medium:
            log.hline()

        self._cache.dump('populations', self.system.numbers - charges)
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
        #return

        funcs = []
        for i in xrange(self._system.natom):
            center = self._system.coordinates[i]
            splines = []
            number = self._system.numbers[i]
            padb_charges = self._proatomdb.get_charges(number, safe=True)
            rtf = self._proatomdb.get_rtransform(number)
            splines = []
            for j0 in xrange(len(padb_charges)-1):
                rad0 = self._proatomdb.get_record(number, padb_charges[j0+1]).rho - \
                       self._proatomdb.get_record(number, padb_charges[j0]).rho
                splines.append(CubicSpline(rad0, rtf=rtf))
                for j1 in xrange(j0+1):
                    rad1 = self._proatomdb.get_record(number, padb_charges[j1+1]).rho - \
                           self._proatomdb.get_record(number, padb_charges[j1]).rho
                    splines.append(CubicSpline(rad0*rad1, rtf=rtf))
            funcs.append((center, splines))
        wcor_fit = self._ui_grid.compute_weight_corrections(funcs, rcond=1.0)
        self._cache.dump('wcor_fit', wcor_fit)

    def _get_constant_fn(self, i, output):
        key = ('isolated_atom', i, 0)
        if key in self._store:
            self._store.load(output, *key)
        else:
            number = self._system.numbers[i]
            center = self._system.coordinates[i]
            spline = self._proatomdb.get_spline(number)
            if log.do_medium:
                log('Computing constant fn %i (n=%i)' % (i, number))
            output[:] = 0
            self._ui_grid.eval_spline(spline, center, output)
            self._store.dump(output, *key)

    def _get_basis_fn(self, i, j, output):
        # A local cache is used that only exists during _init_at_weights:
        key = ('basis', i, j)
        if key in self._store:
            self._store.load(output, *key)
        else:
            number = self._system.numbers[i]
            padb_charges = self._proatomdb.get_charges(number, safe=True)
            center = self._system.coordinates[i]

            if self._proatomdb.get_record(number, padb_charges[0]).pseudo_population == 1:
                if j == 0:
                    spline = self._proatomdb.get_spline(number, {padb_charges[j]: 1})
                    label_q = '%+i' % padb_charges[j]
                else:
                    spline = self._proatomdb.get_spline(number, {padb_charges[j]: 1, padb_charges[j-1]: -1})
                    label_q = '%+i_%+i' % (padb_charges[j], padb_charges[j-1])
            else:
                spline = self._proatomdb.get_spline(number, {padb_charges[j+1]: 1, padb_charges[j]: -1})
                label_q = '%+i_%+i' % (padb_charges[j+1], padb_charges[j])
            if log.do_medium:
                log('Computing basis fn %i_%i (n=%i, q: %s)' % (i, j, number, label_q))
            output[:] = 0
            self._ui_grid.eval_spline(spline, center, output)
            self._store.dump(output, *key)

    def _get_nbasis(self, i):
        number = self._system.numbers[i]
        padb_charges = self._proatomdb.get_charges(number, safe=True)
        result = len(padb_charges)-1
        if self._proatomdb.get_record(number, padb_charges[0]).pseudo_population == 1:
            result += 1
        return result

    def compute_proatom(self, i, output):
        # Get the coefficients for the pro-atom
        nbasis = self._get_nbasis(i)
        pro_coeffs, new = self._cache.load('pro_coeffs', i, alloc=nbasis)

        # Construct the pro-atom
        work = self.zeros()
        self._get_constant_fn(i, output)
        for j in xrange(nbasis):
            if pro_coeffs[j] != 0:
                work[:] = 0
                self._get_basis_fn(i, j, work)
                work *= pro_coeffs[j]
                output += work
        output += 1e-100

    def _update_proatom_pars(self):
        aimdens = self.zeros()
        work0 = self.zeros()
        work1 = self.zeros()
        wcor_fit = self._cache.load('wcor_fit', default=None)
        charges = self._cache.load('charges')
        for i in xrange(self._system.natom):
            # Construct the AIM density
            present = self._store.load(aimdens, 'at_weights', i)
            if not present:
                # construct atomic weight function
                self.compute_proatom(self, i, aimdens)
                aimdens /= self._cache.load('promoldens')
            aimdens *= self._cache.load('moldens')
            self._get_constant_fn(i, work0)
            aimdens -= work0

            # 2) setup equations
            nbasis = self._get_nbasis(i)

            # Preliminary check
            if charges[i] > nbasis:
                raise RuntimeError('The charge on atom %i becomes too positive: %f > %i. (infeasible)' % (i, charges[i], nbasis))

            #    compute A
            A, new = self._cache.load('A', i, alloc=(nbasis, nbasis))
            if new:
                for j0 in xrange(nbasis):
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
            B = np.zeros(nbasis, float)
            for j0 in xrange(nbasis):
                self._get_basis_fn(i, j0, work0)
                B[j0] = self._ui_grid.integrate(aimdens, work0, wcor_fit)
            B /= scales

            C = self._ui_grid.integrate(aimdens, aimdens, wcor_fit)

            # 3) find solution
            #    constraint for total population of pro-atom
            lc_pop = (np.ones(nbasis)/scales, -charges[i])
            #    inequality constraints to keep coefficients larger than -1.
            lcs_par = []
            for j0 in xrange(nbasis):
                lc = np.zeros(nbasis)
                lc[j0] = 1.0/scales[j0]
                lcs_par.append((lc, -1))
                #lcs_par.append((-lc, -1))
            #pro_coeffs = quadratic_solver(A, B, [lc_pop], lcs_par, rcond=0)
            pro_coeffs = quadratic_solver(A, B, [], lcs_par, rcond=0)
            rrmsd = np.sqrt(np.dot(np.dot(A, pro_coeffs) - 2*B, pro_coeffs)/C + 1)

            #    correct for scales
            pro_coeffs /= scales

            if log.do_medium:
                log('            %10i (%.0f%%):&%s' % (i, rrmsd*100, ' '.join('% 6.3f' % c for c in pro_coeffs)))

            self._cache.dump('pro_coeffs', i, pro_coeffs)
