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

from horton.cache import just_once
from horton.grid.int1d import SimpsonIntegrator1D
from horton.grid.cext import CubicSpline, dot_multi
from horton.log import log
from horton.part.hirshfeld import HirshfeldWPart, HirshfeldCPart
from horton.part.hirshfeld_i import HirshfeldIWPart, HirshfeldICPart
from horton.part.linalg import solve_positive, quadratic_solver


__all__ = ['HEBasis', 'HirshfeldEWPart', 'HirshfeldECPart']


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
            rgrid = proatomdb.get_rgrid(number)
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
                    rmsd = np.sqrt(rgrid.integrate(delta, delta))
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
            log('Atom   Z   k label')
            log.hline()
            for i in xrange(len(numbers)):
                for j in xrange(self.get_atom_nbasis(i)):
                    label = self.get_basis_label(i, j)
                    log('%4i %3i %3i %s' % (i, numbers[i], j, label))
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

    def get_constant_lico(self, i):
        return {0: 1}

    def get_basis_rho(self, i, j):
        licos = self.basis_specs[i][2]
        return self.proatomdb.get_rho(self.numbers[i], licos[j])

    def get_basis_spline(self, i, j):
        licos = self.basis_specs[i][2]
        return self.proatomdb.get_spline(self.numbers[i], licos[j])

    def get_basis_lico(self, i, j):
        return self.basis_specs[i][2][j]

    def get_lower_bound(self, i, j):
        lico = self.basis_specs[i][2][j]
        for charge in lico.iterkeys():
            if charge < 0:
                return 0
        return -1

    def get_basis_label(self, i, j):
        licos = self.basis_specs[i][2]
        charges = tuple(sorted(licos[j].keys()))
        if len(charges) == 1:
            return '%+i' % charges
        else:
            return '%+i_%+i' % charges

    def get_initial_propars(self):
        return np.zeros(self.nbasis)

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
    name = 'he'

    def __init__(self, hebasis):
        self._hebasis = hebasis

    def _init_log_scheme(self):
        if log.do_medium:
            log.deflist([
                ('Scheme', 'Hirshfeld-E'),
                ('Convergence threshold', '%.1e' % self._threshold),
                ('Maximum iterations', self._maxiter),
                ('Proatomic DB',  self._proatomdb),
            ])
            log.cite('verstraelen2013', 'the use of Hirshfeld-E partitioning')

    def _get_hebasis(self):
        return self._hebasis

    hebasis = property(_get_hebasis)

    def get_memory_estimates(self):
        if self._greedy:
            return [
                ('Constant', np.ones(self.natom), 0),
                ('Basis', np.array([self.hebasis.get_atom_nbasis(i) for i in xrange(self.natom)]), 0),
            ]
        else:
            return []

    def get_proatom_rho(self, index, propars=None):
        if propars is None:
            propars = self._cache.load('propars')
        begin = self.hebasis.get_atom_begin(index)
        nbasis =  self.hebasis.get_atom_nbasis(index)

        total_lico = {0: 1}
        for j in xrange(nbasis):
            coeff = propars[j+begin]
            lico = self.hebasis.get_basis_lico(index, j)
            for icharge, factor in lico.iteritems():
                total_lico[icharge] = total_lico.get(icharge, 0) + coeff*factor

        number = self._system.numbers[index]
        return self._proatomdb.get_rho(number, total_lico)

    def get_constant(self, index, grid=None):
        spline = self.hebasis.get_constant_spline(index)
        return self.get_somefn(index, spline, ('constant',), 'constant', grid)

    def get_basis(self, index, j, grid=None):
        spline = self.hebasis.get_basis_spline(index, j)
        label = self.hebasis.get_basis_label(index, j)
        return self.get_somefn(index, spline, ('basis', j), 'basis %s' % label, grid)

    def eval_proatom(self, index, output, grid=None):
        # Greedy version of eval_proatom
        output[:] = self.get_constant(index, grid)

        #
        propars = self._cache.load('propars')
        begin = self.hebasis.get_atom_begin(index)
        nbasis =  self.hebasis.get_atom_nbasis(index)

        for j in xrange(nbasis):
            coeff = propars[j+begin]
            if coeff != 0.0:
                output += coeff*self.get_basis(index, j, grid)

        # TODO: check for negative values
        output += 1e-100

    def _init_propars(self):
        self.history_propars = []
        self.history_charges = []
        self.cache.load('charges', alloc=self._system.natom)[0]
        propar_map, propar_names = self.hebasis.get_basis_info()
        self._cache.dump('propar_map', propar_map)
        self._cache.dump('propar_names', np.array(propar_names))
        nbasis = self.hebasis.get_nbasis()
        propars = self.hebasis.get_initial_propars()
        self._cache.dump('propars', propars)
        return propars

    def _update_propars_atom(self, index):
        # Prepare some things
        charges = self._cache.load('charges', alloc=self.system.natom)[0]
        begin = self.hebasis.get_atom_begin(index)
        nbasis = self.hebasis.get_atom_nbasis(index)

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
            lcs_par.append((lc, self.hebasis.get_lower_bound(index, j0)))
        atom_propars = quadratic_solver(A, B, [lc_pop], lcs_par, rcond=0)
        rrms = np.dot(np.dot(A, atom_propars) - 2*B, atom_propars)/C + 1
        if rrms > 0:
            rrmsd = np.sqrt(rrms)
        else:
            rrmsd = -0.01

        #    correct for scales
        atom_propars /= scales

        if log.do_medium:
            log('            %10i (%.0f%%):&%s' % (index, rrmsd*100, ' '.join('% 6.3f' % c for c in atom_propars)))

        self.cache.load('propars')[begin:begin+nbasis] = atom_propars

    def _get_charge_and_delta_aim(self, index):
        # evaluate the constant function, to subtract it from the AIM
        constant = self.get_constant(index)

        # compute delta_aim
        delta_aim = self.get_moldens(index)*self.cache.load('at_weights', index) - constant

        # Integrate out
        grid = self.get_grid(index)
        wcor = self.get_wcor(index)
        charge = -grid.integrate(delta_aim, wcor)
        return charge, delta_aim

    def _get_he_system(self, index, delta_aim):
        number = self.system.numbers[index]
        nbasis = self.hebasis.get_atom_nbasis(index)
        grid = self.get_grid(index)
        wcor_fit = self.get_wcor_fit(index)

        #    Matrix A
        if ('A', number) in self.cache:
            A = self.cache.load('A', number)
        else:
            # Set up system of linear equations:
            A = np.zeros((nbasis, nbasis), float)
            if self.local:
                # In case of local grids, the integration is carried out on
                # a radial grid for efficiency.
                rgrid = self.proatomdb.get_rgrid(number)
                for j0 in xrange(nbasis):
                    rho0 = self.hebasis.get_basis_rho(index, j0)
                    for j1 in xrange(j0+1):
                        rho1 = self.hebasis.get_basis_rho(index, j1)
                        A[j0, j1] = rgrid.integrate(rho0, rho1)
                        A[j1, j0] = A[j0, j1]
            else:
                # In the case of a global grid, the radial integration is not
                # suitable as it does not account for periodic boundary
                # conditions
                basis0 = grid.zeros()
                basis1 = grid.zeros()
                for j0 in xrange(nbasis):
                    basis0[:] = 0.0
                    spline0 = self.hebasis.get_basis_spline(index, j0)
                    self.eval_spline(index, spline0, basis0, label='basis %i' % j0)
                    for j1 in xrange(j0+1):
                        basis1[:] = 0.0
                        spline1 = self.hebasis.get_basis_spline(index, j1)
                        self.eval_spline(index, spline1, basis1, label='basis %i' % j1)
                        A[j0, j1] = grid.integrate(basis0, basis1, wcor_fit)
                        A[j1, j0] = A[j0, j1]

            if (np.diag(A) < 0).any():
                raise ValueError('The diagonal of A must be positive.')

            self.cache.dump('A', number, A)

        #   Matrix B
        B = np.zeros(nbasis, float)
        for j0 in xrange(nbasis):
            basis = self.get_basis(index, j0)
            B[j0] = grid.integrate(delta_aim, basis, wcor_fit)

        #   Constant C
        C = grid.integrate(delta_aim, delta_aim, wcor_fit)

        return A, B, C

    def get_wcor_fit(self, index):
        raise NotImplementedError


class HirshfeldEWPart(HirshfeldEMixin, HirshfeldIWPart):
    def __init__(self, system, grid, proatomdb, local=True, threshold=1e-6, maxiter=500, greedy=False):
        hebasis = HEBasis(system.numbers, proatomdb)
        HirshfeldEMixin.__init__(self, hebasis)
        HirshfeldIWPart.__init__(self, system, grid, proatomdb, local, threshold, maxiter, greedy)

    def get_wcor_fit(self, index):
        return None

    def get_memory_estimates(self):
        return (
            HirshfeldWPart.get_memory_estimates(self) +
            HirshfeldEMixin.get_memory_estimates(self)
        )

    def eval_proatom(self, index, output, grid=None):
        if self._greedy:
            HirshfeldEMixin.eval_proatom(self, index, output, grid)
        else:
            HirshfeldWPart.eval_proatom(self, index, output, grid)


class HirshfeldECPart(HirshfeldEMixin, HirshfeldICPart):
    def __init__(self, system, grid, local, moldens, proatomdb, wcor_numbers, wcor_rcut_max=2.0, wcor_rcond=0.1, threshold=1e-6, maxiter=500, greedy=False):
        '''
           See CPart base class for the description of the arguments.
        '''
        hebasis = HEBasis(system.numbers, proatomdb)
        HirshfeldEMixin.__init__(self, hebasis)
        HirshfeldICPart.__init__(self, system, grid, local, moldens, proatomdb, wcor_numbers, wcor_rcut_max, wcor_rcond, threshold, maxiter, greedy)

    def get_memory_estimates(self):
        if self.local:
            row = [('Weight corrections (fit)', np.array([n in self._wcor_numbers for n in self.system.numbers]), 0)]
        else:
            row = [('Weight corrections (fit)', np.zeros(self.natom), 1)]
        return (
            HirshfeldCPart.get_memory_estimates(self) +
            HirshfeldEMixin.get_memory_estimates(self) +
            row
        )

    def get_wcor_fit_funcs(self, index):
        number = self._system.numbers[index]
        if number not in self.wcor_numbers:
            return []

        center = self._system.coordinates[index]
        atom_nbasis = self.hebasis.get_atom_nbasis(index)
        rtf = self._proatomdb.get_rgrid(self._system.numbers[index]).rtransform
        splines = []
        for j0 in xrange(atom_nbasis):
            rho0 = self.hebasis.get_basis_rho(index, j0)
            splines.append(CubicSpline(rho0, rtf=rtf))
            for j1 in xrange(j0+1):
                rho1 = self.hebasis.get_basis_rho(index, j1)
                splines.append(CubicSpline(rho0*rho1, rtf=rtf))
        return [(center, splines)]

    def get_wcor_fit(self, index=None):
        return self._get_wcor_low('wcor_fit', self.get_wcor_fit_funcs, index)

    def eval_proatom(self, index, output, grid=None):
        if self._greedy:
            HirshfeldEMixin.eval_proatom(self, index, output, grid)
        else:
            HirshfeldCPart.eval_proatom(self, index, output, grid)
