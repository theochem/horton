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

from horton.cext import Cell
from horton.grid.cext import dot_multi, eval_spline_cube
from horton.log import log


__all__ = ['UniformIntGrid']


class UniformIntGrid(object):
    def __init__(self, origin, grid_cell, shape, pbc_active=None):
        if grid_cell.nvec != 3:
            raise ValueError('The cell must be a 3D cell.')
        self.origin = origin
        self.grid_cell = grid_cell
        self.shape = shape
        if pbc_active is None:
            self.pbc_active = np.ones(3, int)
        else:
            self.pbc_active = pbc_active

    def get_cell(self):
        rvecs = (self.grid_cell.rvecs*self.shape.reshape(-1,1))
        return Cell(rvecs[self.pbc_active.astype(bool)])

    def eval_spline(self, spline, center, output):
        eval_spline_cube(spline, center, output, self.origin, self.grid_cell, self.shape, self.pbc_active)

    def integrate(self, *args):
        '''Integrate the product of all arguments

           **Arguments:**

           data1, data2, ...
                All arguments must be arrays with the same size as the number
                of grid points. The arrays contain the functions, evaluated
                at the grid points, that must be multiplied and integrated.
        '''
        # This is often convenient for cube grid data:
        args = [arg.ravel() for arg in args]
        # Similar to conventional integration routine:
        return dot_multi(*args)*self.grid_cell.volume

    def compute_weight_corrections(self, funcs, cache=None):
        '''Computes corrections to the integration weights.

           **Arguments:**

           funcs
                A collection of functions that must integrate exactly with the
                corrected weights. The format is as follows. ``funcs`` is a
                list with tuples that contain three items:

                * center: the center for a set of spherically symmetric
                  functions. In pracice, this will always coincide with th
                  position of a nucleus.

                * rcut: a cutoff radius that limits the corrections to a
                  certain distance from the nucleus. Cutoff spheres of
                  different centers should not overlap.

                * Radial functions specified as a list of tuples. Each tuple
                  contains the following three items:

                    * key: the key used to store the evaluated function in the
                      cache if a cache is provided. If no cache is provided,
                      this may be None.

                    * spline: the radial spline for the spherically symmetric
                      function

                    * integral: the exact integral of the spherically symmetric
                      function [=int_0^r 4*pi*x**2*spline(x)].


           **Optional arguments:**

           cache
                A Cache object in which the evaluated splines are stored to
                avoid their recomputation after the weight corrections are
                computed

           **Return value:**

           The return value is a data array that can be provided as an
           additional argument to the ``integrate`` method. This should
           improve the accuracy of the integration for data that is similar
           to a linear combination of the provided sphericall functions.
        '''
        if cache is None:
            result = np.ones(self.shape, float)
            tmp = np.empty(self.shape, float)
        else:
            result, new = cache.load('wcor', alloc=self.shape)
            result[:] = 1.0
            assert new
        volume = self.grid_cell.volume

        # A safety check to detect overlap of the cutoff spheres.
        # TODO: In the long run, this safety check should be done with a correct
        #       implementation of the MIC.
        cell = self.get_cell()
        for i0 in xrange(len(funcs)):
            center0, rcut0 = funcs[i0][:2]
            for i1 in xrange(i0):
                center1, rcut1 = funcs[i1][:2]
                delta = center1 - center0
                cell.mic(delta)
                dist = np.linalg.norm(delta)
                if dist < rcut0+rcut1:
                    raise ValueError('The cutoff spheres overlap.')

        # Make sure that the rcut spheres fit within the box
        threshold = 0.5*cell.rspacings.min()
        for i0 in xrange(len(funcs)):
            center0, rcut0 = funcs[i0][:2]
            if rcut0 >= threshold:
                raise ValueError('The cutoff spheres have to fit inside the box.')


        icenter = 0
        for center, rcut, rad_funcs in funcs:
            # A) Determine the points inside the cutoff sphere.
            ranges_begin, ranges_end = self.grid_cell.get_ranges_rcut(self.origin, center, rcut)

            # B) Construct a set of grid indexes that lie inside the sphere.
            nselect_max = np.product(ranges_end-ranges_begin)
            indexes = np.zeros((nselect_max, 3), long)
            nselect = self.grid_cell.select_inside(self.origin, center, rcut, ranges_begin, ranges_end, self.shape, self.pbc_active, indexes)
            indexes = indexes[:nselect]

            # C) Allocate the arrays for the least-squares fit of the
            # corrections.
            neq = len(rad_funcs)
            dm = np.zeros((neq+1, nselect), float)
            ev = np.zeros(neq+1, float)

            # D) Fill in the coefficients. This is the expensive part.
            irow = 0
            for key, spline, int_exact in rad_funcs:
                if cache is not None:
                    # If a cache is given, the expensive evaluations are stored
                    # for reuse.
                    tmp, new = cache.load(*key, alloc=self.shape)
                if cache is None or not new:
                    tmp[:] = 0.0
                if log.do_medium:
                    log("Computing spherical function. icenter=%i irow=%i" % (icenter, irow))
                self.eval_spline(spline, center, tmp)
                int_approx = self.integrate(tmp)

                dm[irow] = volume*tmp[indexes[:,0], indexes[:,1], indexes[:,2]]
                ev[irow] = int_exact - int_approx
                irow += 1

            # Also integrate the constant function correctly
            dm[irow] = volume
            ev[irow] = 0.0
            irow += 1

            # rescale equations to optimize condition number
            scales = np.sqrt((dm**2).mean(axis=1))
            dm /= scales.reshape(-1,1)
            ev /= scales

            # E) Find the least norm solution. This part may become more
            # advanced in future.
            U, S, Vt = np.linalg.svd(dm, full_matrices=False)
            assert S[0]*1e-6 < S[-1] # lousy safety check
            corrections = np.dot(Vt.T, np.dot(U.T, ev)/S)

            if log.do_medium:
                rmsd = np.sqrt((corrections**2).mean())
                log('icenter=%i NSELECT=%i CN=%.3e RMSD=%.3e' % (icenter, nselect, S[0]/S[-1], rmsd))

            # F) Fill the corrections into the right place:
            result[indexes[:,0], indexes[:,1], indexes[:,2]] += corrections

            icenter += 1

        if cache is not None:
            cache.discard('tmp')
        return result

    def compute_weight_corrections_brute(self, funcs, cache=None):
        '''Computes corrections to the integration weights.

           **Arguments:**

           funcs
                A collection of functions that must integrate exactly with the
                corrected weights. The format is as follows. ``funcs`` is a
                list with tuples that contain three items:

                * center: the center for a set of spherically symmetric
                  functions. In pracice, this will always coincide with th
                  position of a nucleus.

                * key: the key used to store the evaluated function in the
                  cache if a cache is provided. If no cache is provided,
                  this may be None.

                * spline: the radial spline for the spherically symmetric
                  function

                * integral: the exact integral of the spherically symmetric
                  function [=int_0^r 4*pi*x**2*spline(x)].


           **Optional arguments:**

           cache
                A Cache object in which the evaluated splines are stored to
                avoid their recomputation after the weight corrections are
                computed

           **Return value:**

           The return value is a data array that can be provided as an
           additional argument to the ``integrate`` method. This should
           improve the accuracy of the integration for data that is similar
           to a linear combination of the provided sphericall functions.
        '''
        if cache is None:
            tmp = np.empty(self.shape, float)

        volume = self.grid_cell.volume
        neq = len(funcs)

        # A) Allocate the arrays for the least-squares fit of the
        # corrections.
        npoint = np.product(self.shape)
        dm = np.zeros((neq+1, npoint), float)
        ev = np.zeros(neq+1, float)

        # B) Fill in the rows of the least-squares problem
        for ieq in xrange(neq):
            center, key, spline, int_exact = funcs[ieq]

            if cache is not None:
                # If a cache is given, the expensive evaluations are stored
                # for reuse.
                tmp, new = cache.load(*key, alloc=self.shape)
            if cache is None or not new:
                tmp[:] = 0.0
            if log.do_medium:
                log("Computing spherical function. ieq=%i" % ieq)
            self.eval_spline(spline, center, tmp)
            int_approx = self.integrate(tmp)

            dm[ieq] = volume*tmp.ravel()
            ev[ieq] = int_exact - int_approx

        # Also integrate the constant function correctly
        dm[neq] = 1.0
        ev[neq] = 0.0

        # rescale equations to optimize condition number
        scales = np.sqrt((dm**2).mean(axis=1))
        dm /= scales.reshape(-1,1)
        ev /= scales

        # E) Find the least norm solution. This part may become more
        # advanced in future.
        U, S, Vt = np.linalg.svd(dm, full_matrices=False)
        assert S[0]*1e-6 < S[-1] # lousy safety check
        result = np.dot(Vt.T, np.dot(U.T, ev)/S)

        if log.do_medium:
            rmsd = np.sqrt((result**2).mean())
            log('CN=%.3e RMSD=%.3e' % (S[0]/S[-1], rmsd))

        result += 1

        if cache is not None:
            cache.discard('tmp')
            cache.dump('wcor', result)
        return result
