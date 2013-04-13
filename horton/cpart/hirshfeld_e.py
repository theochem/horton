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

from horton.cpart.hirshfeld_i import HirshfeldICPart
from horton.log import log
from horton.grid.cext import CubicSpline
from horton.dpart.linalg import quadratic_solver


__all__ = [
    'HEBasis', 'HirshfeldECPart',
]


# TODO: reduce duplicate code

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
                self.compute_proatom(i, aimdens)
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
