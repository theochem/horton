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
    'StockholderCPart', 'HirshfeldCPart', 'HirshfeldICPart', 'HirshfeldECPart',
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

    def compute_proatom(self, index, output, window=None):
        # This routine must compute the pro-atom, not from scratch, but based
        # on the info available after the partitioning. The default behavior
        # is to load the pro-atom from the store, but in case the store is fake,
        # this implementation computes the pro-atom on the fly.
        if log.do_medium:
            if window is None:
                shape = self._ui_grid.shape
            else:
                shape = window.shape
            log('Computing proatom %i (%i) with generic routine [%i %i %i]' % (index, self.system.natom, shape[0], shape[1], shape[2]))
        # Construct the pro-atom
        center = self._system.coordinates[index]
        spline = self.get_proatom_spline(index)
        output[:] = 0.0
        if window is None:
            self._ui_grid.eval_spline(spline, center, output)
        else:
            window.eval_spline(spline, center, output)
        output += 1e-100

    def compute_at_weights(self, index, output, window=None):
        self.compute_proatom(index, output, window)
        if window is None:
            output /= self._cache.load('promoldens')
        else:
            promoldens = window.zeros()
            window.extend(self._cache.load('promoldens'), promoldens)
            output /= promoldens

    def get_proatom_spline(self, index):
        raise NotImplementedError


class HirshfeldCPart(StockholderCPart):
    name = 'h'

    def __init__(self, system, ui_grid, moldens, proatomdb, store, smooth=False):
        '''
           See CPart base class for the description of the arguments.
        '''
        self._proatomdb = proatomdb
        CPart.__init__(self, system, ui_grid, moldens, store, smooth)

        # Check of the same (or similar) psuedo potentials were used for the
        # system and the proatoms.
        for i in xrange(self._system.natom):
            number = self._system.numbers[i]
            pseudo_number = self._system.pseudo_numbers[i]
            pseudo_number_expected = self._proatomdb.get_record(number, 0).pseudo_number
            if pseudo_number_expected != pseudo_number:
                raise ValueError('The pseudo number of atom %i does not match with the proatom database (%i!=%i)' % (
                    i, pseudo_number, pseudo_number_expected))

    def _get_proatomdb(self):
        return self._proatomdb

    proatomdb = property(_get_proatomdb)

    def _init_weight_corrections(self):
        funcs = []
        for i in xrange(self._system.natom):
            number = self._system.numbers[i]
            funcs.append((
                self._system.coordinates[i],
                [self._proatomdb.get_spline(number)],
            ))
        wcor = self._ui_grid.compute_weight_corrections(funcs)
        self._cache.dump('wcor', wcor)

    def _init_partitioning(self):
        work = self._ui_grid.zeros()
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

    def get_cutoff_radius(self, index):
        '''The radius at which the weight function goes to zero'''
        rtf = self._proatomdb.get_rtransform(self._system.numbers[index])
        return rtf.radius(rtf.npoint-1)

    def get_proatom_spline(self, index):
        return self._proatomdb.get_spline(self._system.numbers[index])

    def do_all(self):
        names = StockholderCPart.do_all(self)
        self.do_dispersion()
        return names + ['volumes', 'volume_ratios', 'c6s']

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

        volumes, new_volumes = self._cache.load('volumes', alloc=self.system.natom)
        volume_ratios, new_volume_ratios = self._cache.load('volume_ratios', alloc=self.system.natom)
        c6s, new_c6s = self._cache.load('c6s', alloc=self.system.natom)

        if new_volumes or new_volume_ratios or new_c6s:
            self.do_populations()
            self.do_moments()
            k3 = 2
            radial_powers = self._cache.load('radial_powers')
            assert radial_powers[k3] == 3
            radial_moments = self._cache.load('radial_moments')
            populations = self._cache.load('populations')

            if log.do_medium:
                log('Computing atomic dispersion coefficients.')

            for i in xrange(self.system.natom):
                n = self.system.numbers[i]
                if n not in ref_c6s:
                    raise NotImplementedError('No reference C6 value available for atom number %i.' % n)
                volumes[i] = radial_moments[i,k3]/populations[i]
                ref_volume = self.proatomdb.get_record(n, 0).get_moment(3)/n
                volume_ratios[i] = volumes[i]/ref_volume
                c6s[i] = (volume_ratios[i])**2*ref_c6s[n]


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

    def _init_propars(self):
        self.history = []
        return self._cache.load('charges', alloc=self._system.natom)[0]

    def _update_propars(self, charges, work):
        self.history.append(charges.copy())
        for index in xrange(self._system.natom):
            pseudo_population = self.compute_pseudo_population(index, work)
            charges[index] = self.system.pseudo_numbers[index] - pseudo_population

    def _finalize_propars(self, charges):
        self.history.append(charges.copy())
        self._cache.dump('history_charges', np.array(self.history))
        self._cache.dump('populations', self.system.numbers - charges)

    def _init_partitioning(self):
        propars = self._init_propars()
        if log.medium:
            log.hline()
            log('Iteration       Change')
            log.hline()

        counter = 0
        change = 1e100
        work = self._ui_grid.zeros()

        while True:
            counter += 1

            # Update the pro-molecule density
            self._update_promolecule(work)

            # Compute the atomic weight functions if this is useful. This is merely
            # a matter of efficiency.
            self._store_at_weights(work)

            # Update the parameters that determine the pro-atoms.
            old_propars = propars.copy()
            self._update_propars(propars, work)

            # Check for convergence
            change = abs(propars - old_propars).max()
            if log.medium:
                log('%9i   %10.5e' % (counter, change))
            if change < self._threshold or counter >= self._max_iter:
                break

        if log.medium:
            log.hline()

        self._finalize_propars(propars)
        self._cache.dump('niter', counter)
        self._cache.dump('change', change)

    def get_interpolation_info(self, i):
        target_charge = self._cache.load('charges')[i]
        icharge = int(np.floor(target_charge))
        x = target_charge - icharge
        return icharge, x

    def compute_proatom(self, i, output, window=None):
        if self._store.fake or window is not None:
            HirshfeldCPart.compute_proatom(self, i, output, window)
        else:
            # Construct the pro-atom
            icharge, x = self.get_interpolation_info(i)
            self._get_isolated_atom(i, icharge, output)
            if x != 1:
                output *= 1-x
                ipseudo_pop = self.system.pseudo_numbers[i] - icharge
                if ipseudo_pop > 1:
                    tmp = self._ui_grid.zeros()
                    self._get_isolated_atom(i, icharge+1, tmp)
                    tmp *= x
                    output += tmp
            output += 1e-100 # avoid division by zero

    def get_proatom_spline(self, index):
        icharge, x = self.get_interpolation_info(index)
        number = self._system.numbers[index]
        return self._proatomdb.get_spline(number, {icharge: 1-x, icharge+1: x})

    def do_all(self):
        names = HirshfeldCPart.do_all(self)
        return names + ['niter', 'change', 'history_charges']


class HEBasis(object):
    def __init__(self, numbers, proatomdb):
        self.numbers = numbers
        self.proatomdb = proatomdb

        self.nbasis = 0
        self.basis_specs = []

        for i in xrange(len(numbers)):
            number = numbers[i]
            padb_charges = proatomdb.get_charges(number, safe=True)
            complete = proatomdb.get_record(number, padb_charges[0]).pseudo_population == 1
            atom_nbasis = len(padb_charges) - 1 + complete

            licos = []
            self.basis_specs.append([self.nbasis, atom_nbasis, licos])
            for j in xrange(atom_nbasis):
                if complete:
                    if j == 0:
                        licos.append({padb_charges[j]: 1})
                    else:
                        licos.append({padb_charges[j]: 1, padb_charges[j-1]: -1})
                else:
                    licos.append({padb_charges[j+1]: 1, padb_charges[j]: -1})

            self.nbasis += atom_nbasis

    def get_nbasis(self):
        return self.nbasis

    def get_atom_begin(self, i):
        return self.basis_specs[i][0]

    def get_atom_nbasis(self, i):
        return self.basis_specs[i][1]

    def get_constant_spline(self, i):
        return self.proatomdb.get_spline(self.numbers[i])

    def get_basis_spline(self, i, j):
        licos = self.basis_specs[i][2]
        return self.proatomdb.get_spline(self.numbers[i], licos[j])

    def get_basis_lico(self, i, j):
        return self.basis_specs[i][2][j]

    def get_basis_label(self, i, j):
        licos = self.basis_specs[i][2]
        charges = tuple(sorted(licos[j].keys()))
        if len(charges) == 1:
            return '%+i' % charges
        else:
            return '%+i_%+i' % charges

    def get_basis_info(self):
        basis_map = []
        basis_names = []
        for i in xrange(len(self.numbers)):
            begin, nbasis, licos = self.basis_specs[i]
            basis_map.append([begin, nbasis])
            for j in xrange(nbasis):
                basis_names.append('%i:%s' % (i, self.get_basis_label(i, j)))
        return np.array(basis_map), basis_names


class HirshfeldECPart(HirshfeldICPart):
    name = 'he'

    def __init__(self, system, ui_grid, moldens, proatomdb, store, smooth=False, max_iter=100, threshold=1e-4):
        '''
           See CPart base class for the description of the arguments.
        '''
        self._hebasis = HEBasis(system.numbers, proatomdb)
        HirshfeldICPart.__init__(self, system, ui_grid, moldens, proatomdb, store, smooth, max_iter, threshold)


    def _init_weight_corrections(self):
        HirshfeldICPart._init_weight_corrections(self)

        funcs = []
        for i in xrange(self._system.natom):
            center = self._system.coordinates[i]
            splines = []
            atom_nbasis = self._hebasis.get_atom_nbasis(i)
            rtf = self._proatomdb.get_rtransform(self._system.numbers[i])
            splines = []
            for j0 in xrange(atom_nbasis):
                spline0 = self._hebasis.get_basis_spline(i, j0)
                splines.append(spline0)
                for j1 in xrange(j0+1):
                    spline1 = self._hebasis.get_basis_spline(i, j1)
                    splines.append(CubicSpline(spline0.copy_y()*spline1.copy_y(), rtf=rtf))
            funcs.append((center, splines))
        wcor_fit = self._ui_grid.compute_weight_corrections(funcs)
        self._cache.dump('wcor_fit', wcor_fit)

    def _get_constant_fn(self, i, output):
        key = ('isolated_atom', i, 0)
        if key in self._store:
            self._store.load(output, *key)
        else:
            number = self._system.numbers[i]
            center = self._system.coordinates[i]
            spline = self._hebasis.get_constant_spline(i)
            if log.do_medium:
                log('Computing constant fn %i (n=%i)' % (i, number))
            output[:] = 0
            self._ui_grid.eval_spline(spline, center, output)
            self._store.dump(output, *key)

    def _get_basis_fn(self, i, j, output):
        key = ('basis', i, j)
        if key in self._store:
            self._store.load(output, *key)
        else:
            number = self._system.numbers[i]
            padb_charges = self._proatomdb.get_charges(number, safe=True)
            center = self._system.coordinates[i]
            spline = self._hebasis.get_basis_spline(i, j)
            label = self._hebasis.get_basis_label(i, j)
            if log.do_medium:
                log('Computing basis fn %i_%i (n=%i, %s)' % (i, j, number, label))
            output[:] = 0
            self._ui_grid.eval_spline(spline, center, output)
            self._store.dump(output, *key)

    def _init_propars(self):
        self.history = []
        procoeff_map, procoeff_names = self._hebasis.get_basis_info()
        self._cache.dump('procoeff_map', procoeff_map)
        self._cache.dump('procoeff_names', np.array(procoeff_names))
        nbasis = self._hebasis.get_nbasis()
        return self._cache.load('procoeffs', alloc=nbasis)[0]

    def _update_propars(self, procoeffs, work0):
        aimdens = self._ui_grid.zeros()
        work1 = self._ui_grid.zeros()
        wcor = self._cache.load('wcor', default=None)
        wcor_fit = self._cache.load('wcor_fit', default=None)
        charges = self._cache.load('charges', alloc=self._system.natom)[0]
        self.history.append([charges.copy(), procoeffs.copy()])
        for i in xrange(self._system.natom):
            # 1) Construct the AIM density
            present = self._store.load(aimdens, 'at_weights', i)
            if not present:
                # construct atomic weight function
                self.compute_proatom(self, i, aimdens)
                aimdens /= self._cache.load('promoldens')
            aimdens *= self._cache.load('moldens')

            #    compute the charge
            charges[i] = self.system.pseudo_numbers[i] - self._ui_grid.integrate(aimdens, wcor)

            #    subtract the constant function
            self._get_constant_fn(i, work0)
            aimdens -= work0


            # 2) setup equations
            begin = self._hebasis.get_atom_begin(i)
            nbasis = self._hebasis.get_atom_nbasis(i)

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
            atom_procoeffs = quadratic_solver(A, B, [lc_pop], lcs_par, rcond=0)
            rrmsd = np.sqrt(np.dot(np.dot(A, atom_procoeffs) - 2*B, atom_procoeffs)/C + 1)

            #    correct for scales
            atom_procoeffs /= scales

            if log.do_medium:
                log('            %10i (%.0f%%):&%s' % (i, rrmsd*100, ' '.join('% 6.3f' % c for c in atom_procoeffs)))

            procoeffs[begin:begin+nbasis] = atom_procoeffs

    def _finalize_propars(self, procoeffs):
        charges = self._cache.load('charges')
        self.history.append([charges.copy(), procoeffs.copy()])
        self._cache.dump('history_charges', np.array([row[0] for row in self.history]))
        self._cache.dump('history_procoeffs', np.array([row[1] for row in self.history]))
        self._cache.dump('populations', self.system.numbers - charges)

    def compute_proatom(self, i, output, window=None):
        if self._store.fake or window is not None:
            HirshfeldCPart.compute_proatom(self, i, output, window)
        else:
            # Get the coefficients for the pro-atom
            procoeffs = self._cache.load('procoeffs')

            # Construct the pro-atom
            begin = self._hebasis.get_atom_begin(i)
            nbasis =  self._hebasis.get_atom_nbasis(i)
            work = self._ui_grid.zeros()
            self._get_constant_fn(i, output)
            for j in xrange(nbasis):
                if procoeffs[j+begin] != 0:
                    work[:] = 0
                    self._get_basis_fn(i, j, work)
                    work *= procoeffs[j+begin]
                    output += work
            output += 1e-100

    def get_proatom_spline(self, index):
        procoeffs = self._cache.load('procoeffs')
        begin = self._hebasis.get_atom_begin(index)
        nbasis =  self._hebasis.get_atom_nbasis(index)

        total_lico = {0: 1}
        for j in xrange(nbasis):
            coeff = procoeffs[j+begin]
            lico = self._hebasis.get_basis_lico(index, j)
            for icharge, factor in lico.iteritems():
                total_lico[icharge] = total_lico.get(icharge, 0) + coeff*factor

        number = self._system.numbers[index]
        return self._proatomdb.get_spline(number, total_lico)

    def do_all(self):
        names = HirshfeldICPart.do_all(self)
        return names + ['procoeffs', 'procoeff_map', 'procoeff_names', 'history_procoeffs']
