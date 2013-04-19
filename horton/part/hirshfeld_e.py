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
from horton.grid.int1d import SimpsonIntegrator1D
from horton.grid.cext import CubicSpline, dot_multi
from horton.log import log
from horton.part.base import DPart
from horton.part.hirshfeld_i import HirshfeldIDPart, HirshfeldICPart
from horton.part.linalg import solve_positive, quadratic_solver


__all__ = ['HEBasis', 'HirshfeldEDPart', 'HirshfeldECPart']


class HEBasis(object):
    '''Defines the basis set for the promolecule in Hirshfeld-E

       This implementation is based on deviations from the neutral atom. This
       allows one to eliminate basis functions corresponding to very positive
       kations.
    '''
    def __init__(self, numbers, proatomdb):
        self.numbers = numbers
        self.proatomdb = proatomdb

        self.nbasis = 0
        self.basis_specs = []

        for i in xrange(len(numbers)):
            number = numbers[i]
            weights = proatomdb.get_radial_weights(number)
            padb_charges = proatomdb.get_charges(number, safe=True)
            complete = proatomdb.get_record(number, padb_charges[0]).pseudo_population == 1

            licos = []
            atom_nbasis = 0
            for j in xrange(len(padb_charges) - 1 + complete):
                # construct new linear combination
                if complete:
                    if j == 0:
                        lico = {padb_charges[j]: 1}
                    else:
                        lico = {padb_charges[j]: 1, padb_charges[j-1]: -1}
                else:
                    lico = {padb_charges[j+1]: 1, padb_charges[j]: -1}
                # test if this new lico is redundant with respect to previous ones
                redundant = False
                rho = self.proatomdb.get_rho(number, lico)
                for old_lico in licos:
                    old_rho = self.proatomdb.get_rho(number, old_lico)
                    delta = rho - old_rho
                    rmsd = np.sqrt(dot_multi(weights, delta, delta))
                    if rmsd < 1e-8:
                        redundant = True
                        break
                # append if not redundant
                if not redundant:
                    licos.append(lico)
                    atom_nbasis += 1
                else:
                    if log.do_medium:
                        log('Skipping redundant basis function for atom %i: %s' % (i, lico))

            self.basis_specs.append([self.nbasis, atom_nbasis, licos])
            self.nbasis += atom_nbasis

        if log.do_medium:
            log('Hirshfeld-E basis')
            log.hline()
            log('   Z   k label')
            log.hline()
            for i in xrange(len(numbers)):
                for j in xrange(self.get_atom_nbasis(i)):
                    label = self.get_basis_label(i, j)
                    log('%4i %3i %s' % (i, j, label))
            log.hline()

    def get_nbasis(self):
        return self.nbasis

    def get_atom_begin(self, i):
        return self.basis_specs[i][0]

    def get_atom_nbasis(self, i):
        return self.basis_specs[i][1]

    def get_constant_rho(self, i):
        return self.proatomdb.get_rho(self.numbers[i])

    def get_constant_spline(self, i):
        return self.proatomdb.get_spline(self.numbers[i])

    def get_basis_rho(self, i, j):
        licos = self.basis_specs[i][2]
        return self.proatomdb.get_rho(self.numbers[i], licos[j])

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


class HirshfeldEMixin(object):
    def get_proatom_rho(self, index, propars=None):
        if propars is None:
            propars = self._cache.load('propars')
        begin = self._hebasis.get_atom_begin(index)
        nbasis =  self._hebasis.get_atom_nbasis(index)

        total_lico = {0: 1}
        for j in xrange(nbasis):
            coeff = propars[j+begin]
            lico = self._hebasis.get_basis_lico(index, j)
            for icharge, factor in lico.iteritems():
                total_lico[icharge] = total_lico.get(icharge, 0) + coeff*factor

        number = self._system.numbers[index]
        return self._proatomdb.get_rho(number, total_lico)

    def _init_propars(self):
        self.history_propars = []
        self.history_charges = []
        self.cache.load('charges', alloc=self._system.natom)[0]
        propar_map, propar_names = self._hebasis.get_basis_info()
        self._cache.dump('propar_map', propar_map)
        self._cache.dump('propar_names', np.array(propar_names))
        nbasis = self._hebasis.get_nbasis()
        return self._cache.load('propars', alloc=nbasis)[0]

    def _update_propars_atom(self, index):
        # Prepare some things
        charges = self._cache.load('charges', alloc=self.system.natom)[0]
        begin = self._hebasis.get_atom_begin(index)
        nbasis = self._hebasis.get_atom_nbasis(index)

        # Compute charge and delta aim density
        charge, delta_aim = self._get_charge_and_delta_aim(index)
        charges[index] = charge

        # Preliminary check
        if charges[index] > nbasis:
            raise RuntimeError('The charge on atom %i becomes too positive: %f > %i. (infeasible)' % (index, charges[index], nbasis))

        # Define the least-squares system
        A, B, C = self._get_he_system(index, delta_aim)

        # preconditioning
        scales = np.sqrt(np.diag(A))
        A = A/scales/scales.reshape(-1,1)
        B /= scales

        # Find solution
        #    constraint for total population of pro-atom
        lc_pop = (np.ones(nbasis)/scales, -charges[index])
        #    inequality constraints to keep coefficients larger than -1.
        lcs_par = []
        for j0 in xrange(nbasis):
            lc = np.zeros(nbasis)
            lc[j0] = 1.0/scales[j0]
            lcs_par.append((lc, -1))
        atom_propars = quadratic_solver(A, B, [lc_pop], lcs_par, rcond=0)
        rrmsd = np.sqrt(np.dot(np.dot(A, atom_propars) - 2*B, atom_propars)/C + 1)

        #    correct for scales
        atom_propars /= scales

        if log.do_medium:
            log('            %10i (%.0f%%):&%s' % (index, rrmsd*100, ' '.join('% 6.3f' % c for c in atom_propars)))

        self.cache.load('propars')[begin:begin+nbasis] = atom_propars

    def _get_charge_and_delta_aim(self, index):
        raise NotImplementedError

    def _get_he_system(self, index, delta_aim):
        raise NotImplementedError


class HirshfeldEDPart(HirshfeldEMixin, HirshfeldIDPart):
    '''Extended Hirshfeld partitioning'''

    name = 'he'
    options = ['local', 'threshold', 'maxiter']

    def __init__(self, system, grid, proatomdb, local=True, threshold=1e-6, maxiter=500):
        self._hebasis = HEBasis(system.numbers, proatomdb)
        HirshfeldIDPart.__init__(self, system, grid, proatomdb, local, threshold, maxiter)

    def _get_charge_and_delta_aim(self, index):
        # Compute weight
        grid = self.get_grid(index)
        at_weights, new = self.cache.load('at_weights', index, alloc=grid.size)
        self.compute_at_weights(index, at_weights)

        # Compute charge and delta_aim
        work = grid.zeros()
        spline = self._hebasis.get_constant_spline(index)
        center = self.system.coordinates[index]
        grid.eval_spline(spline, center, work)
        delta_aim = self.cache.load('moldens', index)*at_weights - work
        charge = -grid.integrate(delta_aim)
        return charge, delta_aim

    def _get_he_system(self, index, delta_aim):
        number = self.system.numbers[index]
        center = self._system.coordinates[index]
        weights = self.proatomdb.get_radial_weights(number)
        nbasis = self._hebasis.get_atom_nbasis(index)
        grid = self.get_grid(index)

        #    Matrix A
        if self.cache.has('A', number):
            A = self.cache.load('A', number)
        else:
            # Set up system of linear equations:
            A = np.zeros((nbasis, nbasis), float)
            for j0 in xrange(nbasis):
                rho0 = self._hebasis.get_basis_rho(index, j0)
                for j1 in xrange(j0+1):
                    rho1 = self._hebasis.get_basis_rho(index, j1)
                    A[j0, j1] = dot_multi(weights, rho0, rho1)
                    A[j1, j0] = A[j0, j1]

            if (np.diag(A) < 0).any():
                raise ValueError('The diagonal of A must be positive.')

            self.cache.dump('A', number, A)

        #   Matrix B
        B = np.zeros(nbasis, float)
        basis = grid.zeros()
        for j0 in xrange(nbasis):
            basis[:] = 0.0
            spline = self._hebasis.get_basis_spline(index, j0)
            grid.eval_spline(spline, center, basis)
            B[j0] = grid.integrate(delta_aim, basis)

        #   Constant C
        C = grid.integrate(delta_aim, delta_aim)

        return A, B, C

    def do_all(self):
        names = HirshfeldIDPart.do_all(self)
        return names + ['propars', 'propar_map', 'propar_names']


class HirshfeldECPart(HirshfeldEMixin, HirshfeldICPart):
    name = 'he'

    def __init__(self, system, grid, moldens, proatomdb, store, smooth=False, threshold=1e-6, maxiter=500):
        '''
           See CPart base class for the description of the arguments.
        '''
        self._hebasis = HEBasis(system.numbers, proatomdb)
        HirshfeldICPart.__init__(self, system, grid, moldens, proatomdb, store, smooth, threshold, maxiter)

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
                rho0 = self._hebasis.get_basis_rho(i, j0)
                splines.append(CubicSpline(rho0, rtf=rtf))
                for j1 in xrange(j0+1):
                    rho1 = self._hebasis.get_basis_rho(i, j1)
                    splines.append(CubicSpline(rho0*rho1, rtf=rtf))
            funcs.append((center, splines))
        wcor_fit = self.grid.compute_weight_corrections(funcs)
        self._cache.dump('wcor_fit', wcor_fit)

    def _get_constant_fn(self, i, output):
        key = ('isolated_atom', i, 0)
        if key in self._store:
            self._store.load(output, *key)
        else:
            number = self._system.numbers[i]
            spline = self._hebasis.get_constant_spline(i)
            self.compute_spline(i, spline, output, 'n=%i constant' % number)
            self._store.dump(output, *key)

    def _get_basis_fn(self, i, j, output):
        key = ('basis', i, j)
        if key in self._store:
            self._store.load(output, *key)
        else:
            number = self._system.numbers[i]
            spline = self._hebasis.get_basis_spline(i, j)
            label = self._hebasis.get_basis_label(i, j)
            self.compute_spline(i, spline, output, 'n=%i %s' % (number, label))
            self._store.dump(output, *key)

    def _get_charge_and_delta_aim(self, index):
        delta_aim = self.grid.zeros()
        work = self.cache.load('work0', alloc=self.grid.shape)[0]
        wcor = self._cache.load('wcor', default=None)

        # 1) Construct the AIM density
        present = self._store.load(delta_aim, 'at_weights', index)
        if not present:
            # construct atomic weight function
            self.compute_proatom(index, delta_aims)
            delta_aim /= self._cache.load('promoldens')
        delta_aim *= self._cache.load('moldens')

        #    subtract the constant function
        self._get_constant_fn(index, work)
        delta_aim -= work

        # compute the charge
        charge = -self.grid.integrate(delta_aim, wcor)

        return charge, delta_aim

    def _get_he_system(self, index, delta_aim):
        work0 = self.cache.load('work0', alloc=self.grid.shape)[0]
        work1 = self.cache.load('work1', alloc=self.grid.shape)[0]
        wcor_fit = self._cache.load('wcor_fit', default=None)
        nbasis = self._hebasis.get_atom_nbasis(index)

        #    compute A
        A, new = self._cache.load('A', index, alloc=(nbasis, nbasis))
        if new:
            for j0 in xrange(nbasis):
                self._get_basis_fn(index, j0, work0)
                for j1 in xrange(j0+1):
                    self._get_basis_fn(index, j1, work1)
                    A[j0,j1] = self.grid.integrate(work0, work1, wcor_fit)
                    A[j1,j0] = A[j0,j1]

        #    compute B
        B = np.zeros(nbasis, float)
        for j0 in xrange(nbasis):
            self._get_basis_fn(index, j0, work0)
            B[j0] = self.grid.integrate(delta_aim, work0, wcor_fit)

        #    constant c
        C = self.grid.integrate(delta_aim, delta_aim, wcor_fit)

        return A, B, C

    def compute_proatom(self, i, output, window=None):
        if self._store.fake or window is not None:
            HirshfeldICPart.compute_proatom(self, i, output, window)
        else:
            # Get the coefficients for the pro-atom
            propars = self._cache.load('propars')

            # Construct the pro-atom
            begin = self._hebasis.get_atom_begin(i)
            nbasis =  self._hebasis.get_atom_nbasis(i)
            work = self.cache.load('work1', alloc=self.grid.shape)[0]
            self._get_constant_fn(i, output)
            for j in xrange(nbasis):
                if propars[j+begin] != 0:
                    self._get_basis_fn(i, j, work)
                    work *= propars[j+begin]
                    output += work
            output += 1e-100

    def do_all(self):
        names = HirshfeldICPart.do_all(self)
        return names + ['propars', 'propar_map', 'propar_names']
