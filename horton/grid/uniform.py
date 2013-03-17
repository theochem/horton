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
from horton.grid.int1d import SimpsonIntegrator1D
from horton.grid.cext import dot_multi, dot_multi_poly_cube, eval_spline_cube
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
            self.pbc_active = pbc_active.astype(int)

    @classmethod
    def from_hdf5(cls, grp, lf):
        return cls(
            grp['origin'][:],
            Cell.from_hdf5(grp['grid_cell'], lf),
            grp['shape'][:],
            grp['pbc_active'][:],
        )

    def to_hdf5(self, grp):
        subgrp = grp.require_group('grid_cell')
        self.grid_cell.to_hdf5(subgrp)
        grp['origin'] = self.origin
        grp['shape'] = self.shape
        grp['pbc_active'] = self.pbc_active

    def get_cell(self):
        rvecs = (self.grid_cell.rvecs*self.shape.reshape(-1,1))
        return Cell(rvecs[self.pbc_active.astype(bool)])

    def _get_size(self):
        return np.product(self.shape)

    size = property(_get_size)

    def eval_spline(self, spline, center, output):
        eval_spline_cube(spline, center, output, self.origin, self.grid_cell, self.shape, self.pbc_active)

    def integrate(self, *args, **kwargs):
        '''Integrate the product of all arguments

           **Arguments:**

           data1, data2, ...
                All arguments must be arrays with the same size as the number
                of grid points. The arrays contain the functions, evaluated
                at the grid points, that must be multiplied and integrated.

           **Optional arguments:**

           TODO
        '''
        center = kwargs.pop('center', None)
        mask = kwargs.pop('mask', True)
        powx = kwargs.pop('powx', 0)
        powy = kwargs.pop('powy', 0)
        powz = kwargs.pop('powz', 0)
        powr = kwargs.pop('powr', 0)
        if len(kwargs) > 0:
            raise TypeError('Unexpected arguments: %s' % tuple(kwargs.keys()))
        # This is often convenient for cube grid data:
        args = [arg.ravel() for arg in args if arg is not None]
        # Similar to conventional integration routine:
        if center is None:
            return dot_multi(*args)*self.grid_cell.volume
        else:
            return dot_multi_poly_cube(
                args, self.origin, self.grid_cell, self.shape, self.pbc_active,
                self.get_cell(), center, mask, powx, powy, powz, powr)*self.grid_cell.volume

    def compute_weight_corrections(self, funcs, rcut_scale=0.9, rcut_max=2.0, rcond=0.1):
        '''Computes corrections to the integration weights.

           **Arguments:**

           funcs
                A collection of functions that must integrate exactly with the
                corrected weights. The format is as follows. ``funcs`` is a
                list with tuples that contain three items:

                * center: the center for a set of spherically symmetric
                  functions. In pracice, this will always coincide with th
                  position of a nucleus.

                * Radial functions specified as a list of splines.

           **Optional arguments:**

           rcut_scale
                For center (of a spherical function), radii of non-overlapping
                spheres are determined by setting the radius of each sphere at
                0.5*rcut_scale*(distance to nearest atom or periodic image).

           rcut_max
                To avoid gigantic cutoff spheres, one may use rcut_max to set
                the maximum radius of the cutoff sphere.

           rcond
                The regulatization strength for the weight correction equations.
                This should not be too low. Current value is a compromise
                between accuracy and transferability of the weight corrections.

           **Return value:**

           The return value is a data array that can be provided as an
           additional argument to the ``integrate`` method. This should
           improve the accuracy of the integration for data that is similar
           to a linear combination of the provided sphericall functions.
        '''
        result = np.ones(self.shape, float)
        volume = self.grid_cell.volume

        # initialize cutoff radii
        cell = self.get_cell()
        rcut_max = min(rcut_max, 0.5*rcut_scale*cell.rspacings.min())
        rcuts = np.zeros(len(funcs)) + rcut_max

        # determine safe cutoff radii
        for i0 in xrange(len(funcs)):
            center0, rcut0 = funcs[i0][:2]
            for i1 in xrange(i0):
                center1, rcut1 = funcs[i1][:2]
                delta = center1 - center0
                cell.mic(delta)
                dist = np.linalg.norm(delta)
                rcut = 0.5*rcut_scale*dist
                rcuts[i0] = min(rcut, rcuts[i0])
                rcuts[i1] = min(rcut, rcuts[i1])

        def get_aux_grid(center, aux_rcut):
            ranges_begin, ranges_end = self.grid_cell.get_ranges_rcut(center-self.origin, aux_rcut)
            aux_origin = self.origin.copy()
            self.grid_cell.add_vec(aux_origin, ranges_begin)
            aux_shape = ranges_end - ranges_begin
            aux_grid = UniformIntGrid(aux_origin, self.grid_cell, aux_shape, pbc_active=np.zeros(3, float))
            return aux_grid, -ranges_begin

        def get_tapered_spline(spline, rcut, aux_rcut):
            assert rcut < aux_rcut
            rtf = spline.rtransform
            r = rtf.get_radii()
            # Get original spline stuff
            y = spline.copy_y()
            d = spline.deriv(r)
            # adapt cutoffs to indexes of radial grid
            ibegin = r.searchsorted(rcut)
            iend = r.searchsorted(aux_rcut)-1
            rcut = r[ibegin]
            aux_rcut = r[iend]
            # define tapering function
            sy = np.zeros(len(r))
            sd = np.zeros(len(r))
            sy[:ibegin] = 1.0
            scale = np.pi/(aux_rcut-rcut)
            x = scale*(r[ibegin:iend+1]-rcut)
            sy[ibegin:iend+1] = 0.5*(np.cos(x)+1)
            sd[ibegin:iend+1] = -0.5*(np.sin(x))*scale
            # construct product
            ty = y*sy
            td = d*sy+y*sd
            # construct spline
            tapered_spline = CubicSpline(ty, td, rtf)
            # compute integral
            int1d = SimpsonIntegrator1D()
            int_exact = 4*np.pi*dot_multi(
                ty, r, r, int1d.get_weights(len(r)), rtf.get_volume_elements()
            )
            # done
            return tapered_spline, int_exact

        icenter = 0
        for (center, splines), rcut in zip(funcs, rcuts):
            # A) Determine the points inside the cutoff sphere.
            ranges_begin, ranges_end = self.grid_cell.get_ranges_rcut(center-self.origin, rcut)

            # B) Construct a set of grid indexes that lie inside the sphere.
            nselect_max = np.product(ranges_end-ranges_begin)
            indexes = np.zeros((nselect_max, 3), int)
            nselect = self.grid_cell.select_inside(self.origin, center, rcut, ranges_begin, ranges_end, self.shape, self.pbc_active, indexes)
            indexes = indexes[:nselect]

            # C) Set up an integration grid for the tapered spline
            aux_rcut = 2*rcut
            aux_grid, aux_offset = get_aux_grid(center, aux_rcut)
            aux_indexes = (indexes + aux_offset) % self.shape

            # D) Allocate the arrays for the least-squares fit of the
            # corrections.
            neq = len(splines)
            dm = np.zeros((neq+1, nselect), float)
            ev = np.zeros(neq+1, float)

            # E) Fill in the coefficients. This is the expensive part.
            ieq = 0
            tmp = np.zeros(aux_grid.shape)
            av = np.zeros(neq+1)
            for spline in splines:
                if log.do_medium:
                    log("Computing spherical function. icenter=%i ieq=%i" % (icenter, ieq))
                tapered_spline, int_exact = get_tapered_spline(spline, rcut, aux_rcut)
                tmp[:] = 0.0
                aux_grid.eval_spline(tapered_spline, center, tmp)
                int_approx = self.integrate(tmp)
                dm[ieq] = volume*tmp[aux_indexes[:,0], aux_indexes[:,1], aux_indexes[:,2]]
                av[ieq] = int_approx
                ev[ieq] = int_exact - int_approx
                ieq += 1

            # Add error on constant function
            dm[neq] = volume
            av[neq] = 0.0
            ev[neq] = 0.0

            # rescale equations to optimize condition number
            scales = np.sqrt((dm**2).mean(axis=1))
            dm /= scales.reshape(-1,1)
            ev /= scales
            av /= scales

            # E) Find a regularized least norm solution.
            U, S, Vt = np.linalg.svd(dm, full_matrices=False)
            ridge = rcond*S[0]
            Sinv = S/(ridge**2+S**2)
            corrections = np.dot(Vt.T, np.dot(U.T, ev)*Sinv)

            # constrain the solution to integrate constant function exactly
            HtVSinv = Vt.sum(axis=1)*Sinv
            mu = corrections.sum()/np.dot(HtVSinv,HtVSinv)
            corrections -= mu*np.dot(Vt.T, Vt.sum(axis=1)*Sinv**2)

            if log.do_medium:
                rmsd = np.sqrt((corrections**2).mean())
                log('icenter=%i NSELECT=%i CN=%.3e RMSD=%.3e' % (icenter, nselect, S[0]/S[-1], rmsd))
                mv = np.dot(dm, corrections)
                for ieq in xrange(neq):
                    log('   spline %3i    error %+.3e   orig %+.3e  exact %+.3e' % (ieq, ev[ieq]-mv[ieq], ev[ieq], ev[ieq]+av[ieq]))
                log('   constant      error %+.3e' % corrections.sum())

            # F) Fill the corrections into the right place:
            result[indexes[:,0], indexes[:,1], indexes[:,2]] += corrections

            icenter += 1

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
