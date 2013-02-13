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

from horton.cpart.base import CPart, CCPart
from horton.cache import just_once, Cache
from horton.log import log


__all__ = [
    'HirshfeldCPart', 'HirshfeldICPart', 'HirshfeldECPart',
    'HirshfeldCCPart', 'HirshfeldICCPart', 'HirshfeldECCPart',
]


# TODO: maxiter, tolerance parameter


class HirshfeldCPart(CPart):
    def __init__(self, system, ui_grid, mol_dens, proatomdb):
        '''
           **Arguments:**

           system
                The system to be partitioned.

           ui_grid
                The uniform integration grid based on the cube file.

           mol_dens
                The all-electron molecular density grid data.

           proatomdb
                The database of proatomic densities.
        '''
        self._proatomdb = proatomdb
        CPart.__init__(self, system, ui_grid, mol_dens)

    def _get_proatomdb(self):
        return self._proatomdb

    proatomdb = property(_get_proatomdb)

    def _compute_pro_atom(self, i):
        if log.do_medium:
            log('Computing pro-atom %i (%i)' % (i, self._system.natom))
        # Pro-atoms are (temporarily) stored in at_weights for efficiency.
        pro_atom, new = self._cache.load('at_weights', i, alloc=self._ui_grid.shape)
        # Construct the pro-atom
        assert new
        number = self._system.numbers[i]
        center = self._system.coordinates[i]
        spline = self._proatomdb.get_hirshfeld_proatom_fn(number)
        self._ui_grid.eval_spline(spline, center, pro_atom)
        pro_atom += 1e-100 # avoid division by zero

    def _compute_pro_molecule(self):
        # The pro-molecule is (temporarily) stored in rel_mol_dens for efficiency.
        pro_molecule, new = self._cache.load('rel_mol_dens', alloc=self._ui_grid.shape)
        if not new:
            pro_molecule[:] = 0
        for i in xrange(self._system.natom):
            pro_atom = self._cache.load('at_weights', i)
            pro_molecule += pro_atom

    def _compute_at_weights(self):
        pro_molecule = self._cache.load('rel_mol_dens')
        for i in xrange(self._system.natom):
            pro_atom = self._cache.load('at_weights', i)
            pro_atom /= pro_molecule

    def _compute_rel_dens(self):
        rel_mol_dens = self._cache.load('rel_mol_dens')
        rel_mol_dens *= -1
        rel_mol_dens += self._cache.load('mol_dens')

    @just_once
    def _init_at_weights(self):
        # A) Compute all the pro-atomic densities.
        for i in xrange(self._system.natom):
            self._compute_pro_atom(i)
        self._compute_pro_molecule()
        self._compute_at_weights()
        self._compute_rel_dens()
        # Store reference populations
        self._cache.dump('ref_populations', self._system.numbers.astype(float))


class HirshfeldICPart(HirshfeldCPart):
    def _get_isolated_atom(self, i, pop):
        # A local cache is used that only exists during _init_at_weights:
        isolated_atom, new = self._local_cache.load('isolated_atom', i, pop, alloc=self._ui_grid.shape)
        if new:
            number = self._system.numbers[i]
            center = self._system.coordinates[i]
            spline = self._proatomdb.get_hirshfeld_i_proatom_fn(number, pop)
            if log.do_medium:
                log('Computing isolated atom %i (n=%i, pop=%i)' % (i, number, pop))
            self._ui_grid.eval_spline(spline, center, isolated_atom)
        return isolated_atom

    def _compute_pro_atom(self, i, pop):
        # Pro-atoms are (temporarily) stored in at_weights for efficiency.
        pro_atom, new = self._cache.load('at_weights', i, alloc=self._ui_grid.shape)
        # Construct the pro-atom
        ipop = int(np.floor(pop))
        pro_atom[:] = self._get_isolated_atom(i, ipop)
        if float(ipop) != pop:
            x = pop - ipop
            pro_atom *= 1-x
            tmp = self._local_cache.load('tmp', alloc=self._ui_grid.shape)[0]
            tmp[:] = self._get_isolated_atom(i, ipop+1)
            tmp *= x
            pro_atom += tmp
        pro_atom += 1e-100 # avoid division by zero
        return pro_atom

    @just_once
    def _init_at_weights(self):
        self._local_cache = Cache()

        ref_populations = self._system.numbers.astype(float)
        if log.medium:
            log.hline()
            log('Iteration       Change')
            log.hline()
        counter = 0
        while True:
            # Construct pro-atoms
            for i in xrange(self._system.natom):
                self._compute_pro_atom(i, ref_populations[i])
            self._compute_pro_molecule()
            self._compute_at_weights()
            self._compute_rel_dens()
            new_populations = self._compute_rel_populations() + ref_populations
            change = abs(new_populations - ref_populations).max()
            if log.medium:
                log('%9i   %10.5e' % (counter, change))
            if change < 1e-4:
                break

            ref_populations = new_populations
            counter += 1

        if log.medium:
            log.hline()

        self._cache.dump('populations', new_populations)
        self._cache.dump('ref_populations', ref_populations)
        del self._local_cache


def positive_solve(A, B, constraint_pars=None):
    from scipy.optimize import fmin_slsqp
    def cost_fn(x):
        return np.dot(x, np.dot(A, x)) - 2*np.dot(B, x)

    def cost_fn_gradient(x):
        return 2*(np.dot(A, x) - B)

    if constraint_pars is None:
        f_eqcons = None
        fprime_eqcons = None
    else:
        pop, pop_begin, pop_end = constraint_pars
        pops = np.arange(pop_begin, pop_end)

        def f_eqcons(x):
            return np.array([np.dot(x, pops) - pop])

        def fprime_eqcons(x):
            return np.array([pops])

    N = len(B)
    x0 = np.ones(N, float)/N
    x1 = fmin_slsqp(cost_fn, x0, bounds=[(0, 10)]*N, fprime=cost_fn_gradient, f_eqcons=f_eqcons, fprime_eqcons=fprime_eqcons, iprint=0, acc=1e-10)
    return x1


class HirshfeldECPart(HirshfeldICPart):
    def _compute_pro_atom(self, i):
        # Pro-atoms are (temporarily) stored in at_weights for efficiency.
        pro_atom, new = self._cache.load('at_weights', i, alloc=self._ui_grid.shape)
        number = self._system.numbers[i]
        if new:
            # just use the Hirshfeld definition to get started
            pro_atom[:] = self._get_isolated_atom(i, number)
            return number
        else:
            # do a least-squares fit to the previos AIM

            # 1) construct aim in temporary array
            rho_aim = self._local_cache.load('rho_aim', alloc=self._ui_grid.shape)[0]
            rho_aim[:] = pro_atom
            rho_aim *= self._cache.load('mol_dens')

            # 2) setup equations
            pop_begin, pop_end = self._proatomdb.get_pop_range(number)
            neq = pop_end - pop_begin
            B = np.zeros(neq, float)
            for j0 in xrange(neq):
                pop0 = pop_begin + j0
                rho0_pop0 = self._get_isolated_atom(i, pop0)
                B[j0] = self._ui_grid.integrate(rho_aim, rho0_pop0)

            A, new = self._local_cache.load('fit_A', i, alloc=(neq, neq))
            if new:
                for j0 in xrange(neq):
                    pop0 = pop_begin + j0
                    rho0_pop0 = self._get_isolated_atom(i, pop0)
                    for j1 in xrange(j0+1):
                        pop1 = pop_begin + j1
                        rho0_pop1 = self._get_isolated_atom(i, pop1)
                        A[j0,j1] = self._ui_grid.integrate(rho0_pop0, rho0_pop1)
                        A[j1,j0] = A[j0,j1]

            # 3) find positive solution
            c = positive_solve(A, B)

            # 4) construct the pro-atom
            tmp = rho_aim
            del rho_aim
            pro_atom[:] = 0
            ref_population = 0
            for j0 in xrange(neq):
                pop0 = pop_begin + j0
                tmp[:] = self._get_isolated_atom(i, pop0)
                tmp *= c[j0]
                pro_atom += tmp
                ref_population += c[j0]*pop0

            return ref_population

    @just_once
    def _init_at_weights(self):
        self._local_cache = Cache()

        if log.medium:
            log.hline()
            log('Iteration       Change')
            log.hline()
        counter = 0
        ref_populations = np.zeros(self._system.natom)
        old_populations = self._system.numbers.astype(float)
        while True:
            # Construct pro-atoms
            for i in xrange(self._system.natom):
                ref_populations[i] = self._compute_pro_atom(i)
            self._compute_pro_molecule()
            self._compute_at_weights()
            self._compute_rel_dens()
            new_populations = self._compute_rel_populations() + ref_populations
            change = abs(new_populations - old_populations).max()
            if log.medium:
                log('%9i   %10.5e' % (counter, change))
            if change < 1e-4:
                break

            old_populations = new_populations
            counter += 1

        if log.medium:
            log.hline()

        self._cache.dump('populations', new_populations)
        self._cache.dump('ref_populations', ref_populations)
        del self._local_cache


class HirshfeldCCPart(CCPart):
    def __init__(self, system, ui_grid, mol_dens, proatomdb):
        '''
           **Arguments:**

           system
                The system to be partitioned.

           ui_grid
                The uniform integration grid based on the cube file.

           mol_dens
                The all-electron molecular density grid data.

           proatomdb
                The database of proatomic densities.
        '''
        self._proatomdb = proatomdb
        CCPart.__init__(self, system, ui_grid, mol_dens)

    def _get_proatomdb(self):
        return self._proatomdb

    proatomdb = property(_get_proatomdb)

    @just_once
    def _init_weight_corrections(self):
        if True:
            funcs = []
            for i in xrange(self._system.natom):
                n = self._system.numbers[i]
                funcs.append((
                    self._system.coordinates[i],
                    1.0, # TODO determine radii in clever way or write new correction algo
                    [(('isolated_atom', i, n), self._proatomdb.get_hirshfeld_proatom_fn(n), float(n))],
                ))
            self._ui_grid.compute_weight_corrections(funcs, cache=self._cache)
        else:
            # TODO: investigate this option in more detail
            funcs = []
            for i in xrange(self._system.natom):
                n = self._system.numbers[i]
                pop_begin, pop_end = self._proatomdb.get_pop_range(n)
                for pop in xrange(pop_begin, pop_end):
                    funcs.append((
                        self._system.coordinates[i],
                        ('isolated_atom', i, pop),
                        self._proatomdb.get_hirshfeld_i_proatom_fn(n, pop),
                        float(n),
                    ))
            self._ui_grid.compute_weight_corrections_brute(funcs, cache=self._cache)

    def _compute_pro_atom(self, i):
        n = self._system.numbers[i]
        self._cache.rename(('isolated_atom', i, n), ('at_weights', i))
        pro_atom = self._cache.load('at_weights', i)
        pro_atom += 1e-100

    def _compute_pro_molecule(self):
        pro_molecule, new = self._cache.load('pro_molecule', alloc=self._ui_grid.shape)
        if not new:
            pro_molecule[:] = 0
        for i in xrange(self._system.natom):
            pro_atom = self._cache.load('at_weights', i)
            pro_molecule += pro_atom

    def _compute_at_weights(self):
        pro_molecule = self._cache.load('pro_molecule')
        for i in xrange(self._system.natom):
            pro_atom = self._cache.load('at_weights', i)
            pro_atom /= pro_molecule

    @just_once
    def _init_at_weights(self):
        # A) Compute all the pro-atomic densities.
        for i in xrange(self._system.natom):
            self._compute_pro_atom(i)
        self._compute_pro_molecule()
        self._compute_at_weights()
        self._cache.discard('pro_molecule')



class HirshfeldICCPart(HirshfeldCCPart):
    def _get_isolated_atom(self, i, pop):
        # A local cache is used that only exists during _init_at_weights:
        isolated_atom, new = self._cache.load('isolated_atom', i, pop, alloc=self._ui_grid.shape)
        if new:
            number = self._system.numbers[i]
            center = self._system.coordinates[i]
            spline = self._proatomdb.get_hirshfeld_i_proatom_fn(number, pop)
            if log.do_medium:
                log('Computing isolated atom %i (n=%i, pop=%i)' % (i, number, pop))
            self._ui_grid.eval_spline(spline, center, isolated_atom)
        return isolated_atom

    def _compute_pro_atom(self, i, pop):
        # Pro-atoms are (temporarily) stored in at_weights for efficiency.
        pro_atom, new = self._cache.load('at_weights', i, alloc=self._ui_grid.shape)
        # Construct the pro-atom
        ipop = int(np.floor(pop))
        pro_atom[:] = self._get_isolated_atom(i, ipop)
        if float(ipop) != pop:
            x = pop - ipop
            pro_atom *= 1-x
            tmp = self._local_cache.load('tmp', alloc=self._ui_grid.shape)[0]
            tmp[:] = self._get_isolated_atom(i, ipop+1)
            tmp *= x
            pro_atom += tmp
        pro_atom += 1e-100 # avoid division by zero
        return pro_atom

    @just_once
    def _init_at_weights(self):
        self._local_cache = Cache()

        old_populations = self._system.numbers.astype(float)
        if log.medium:
            log.hline()
            log('Iteration       Change')
            log.hline()
        counter = 0
        while True:
            # Construct pro-atoms
            for i in xrange(self._system.natom):
                self._compute_pro_atom(i, old_populations[i])
            self._compute_pro_molecule()
            self._compute_at_weights()
            new_populations = self._compute_populations()
            change = abs(new_populations - old_populations).max()
            if log.medium:
                log('%9i   %10.5e' % (counter, change))
            if change < 1e-4:
                break

            old_populations = new_populations
            counter += 1

        if log.medium:
            log.hline()

        self._cache.dump('populations', new_populations)
        del self._local_cache


class HirshfeldECCPart(HirshfeldICCPart):
    @just_once
    def _init_weight_corrections(self):
        funcs = []
        from horton.grid import SimpsonIntegrator1D, CubicSpline
        int1d = SimpsonIntegrator1D()
        rtf = self._proatomdb._rtransform
        radii = rtf.get_radii()
        weights1d = (4*np.pi) * radii**2 * rtf.get_volume_elements() * int1d.get_weights(len(radii))

        for i in xrange(self._system.natom):
            n = self._system.numbers[i]
            record = self._proatomdb._records[(n, n)]
            spline1 = CubicSpline(record, rtf=rtf)
            spline2 = CubicSpline(record**2, rtf=rtf)
            int_exact1 = float(n)
            int_exact2 = np.dot(record**2, weights1d)
            funcs.append((
                self._system.coordinates[i],
                1.0, # TODO determine radii in clever way or write new correction algo
                [
                    (('isolated_atom', i, n), spline1, int_exact1),
                    (('tmp'), spline2, int_exact2),
                ],
            ))
        self._ui_grid.compute_weight_corrections(funcs, cache=self._cache)

    def _compute_pro_atom(self, i, pop):
        # Pro-atoms are (temporarily) stored in at_weights for efficiency.
        wcor = self._cache.load('wcor')
        pro_atom, new = self._cache.load('at_weights', i, alloc=self._ui_grid.shape)
        number = self._system.numbers[i]
        if new:
            # just use the Hirshfeld definition to get started
            pro_atom[:] = self._get_isolated_atom(i, number)
            return number
        else:
            # do a least-squares fit to the previos AIM

            # 1) construct aim in temporary array
            rho_aim = self._local_cache.load('rho_aim', alloc=self._ui_grid.shape)[0]
            rho_aim[:] = pro_atom
            rho_aim *= self._cache.load('mol_dens')

            # 2) setup equations
            pop_begin, pop_end = self._proatomdb.get_pop_range(number)
            neq = pop_end - pop_begin
            B = np.zeros(neq, float)
            for j0 in xrange(neq):
                pop0 = pop_begin + j0
                rho0_pop0 = self._get_isolated_atom(i, pop0)
                B[j0] = self._ui_grid.integrate(rho_aim, rho0_pop0, wcor)

            A, new = self._local_cache.load('fit_A', i, alloc=(neq, neq))
            if new:
                for j0 in xrange(neq):
                    pop0 = pop_begin + j0
                    rho0_pop0 = self._get_isolated_atom(i, pop0)
                    for j1 in xrange(j0+1):
                        pop1 = pop_begin + j1
                        rho0_pop1 = self._get_isolated_atom(i, pop1)
                        A[j0,j1] = self._ui_grid.integrate(rho0_pop0, rho0_pop1, wcor)
                        A[j1,j0] = A[j0,j1]

            # 3) find positive solution
            constraint_pars = None
            #constraint_pars = (pop, pop_begin, pop_end)
            c = positive_solve(A, B, constraint_pars)

            # 4) construct the pro-atom
            tmp = rho_aim
            del rho_aim
            pro_atom[:] = 0
            for j0 in xrange(neq):
                pop0 = pop_begin + j0
                tmp[:] = self._get_isolated_atom(i, pop0)
                tmp *= c[j0]
                pro_atom += tmp

    @just_once
    def _init_at_weights(self):
        self._local_cache = Cache()

        if log.medium:
            log.hline()
            log('Iteration       Change')
            log.hline()
        counter = 0
        old_populations = self._system.numbers.astype(float)
        while True:
            # Construct pro-atoms
            for i in xrange(self._system.natom):
                self._compute_pro_atom(i, old_populations[i])
            self._compute_pro_molecule()
            self._compute_at_weights()
            new_populations = self._compute_populations()
            change = abs(new_populations - old_populations).max()
            if log.medium:
                log('%9i   %10.5e' % (counter, change))
            if change < 1e-4:
                break

            old_populations = new_populations
            counter += 1

        if log.medium:
            log.hline()

        self._cache.dump('populations', new_populations)
        del self._local_cache
