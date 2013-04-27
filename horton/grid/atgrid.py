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


from horton.context import context
from horton.grid.base import IntGrid
from horton.grid.cext import lebedev_laikov_npoints, lebedev_laikov_sphere, RTransform, dot_multi_parts
from horton.grid.int1d import SimpsonIntegrator1D
from horton.grid.radial import RadialIntGrid
from horton.log import log


__all__ = [
    'AtomicGrid', 'interpret_atspec', 'get_random_rotation', 'ATGridFamily',
    'atgrid_families',
]


class AtomicGrid(IntGrid):
    def __init__(self, number, center, atspec='tv-13.1-3', random_rotate=True, points=None):
        '''
           **Arguments:**

           center
                The center of the radial grid

           atspec
                A specifications of the atomic grid. (See interpret_atspec for
                more details.)

           **Optional arguments:**

           random_rotate
                When set to False, the random rotation of the grid points is
                disabled. Such random rotation improves the accuracy of the
                integration, but leads to small random changes in the results
                that are not reproducible.

           points
                Array to store the grid points
        '''
        self._number = number
        self._center = center
        self._rtransform, self._int1d, self._nlls = interpret_atspec(number, atspec)
        self._random_rotate = random_rotate
        self._rgrid = RadialIntGrid(self._rtransform, self._int1d)

        size = self._nlls.sum()
        if points is None:
            points = np.zeros((size, 3), float)
        else:
            assert len(points) == size
        weights = np.zeros(size, float)
        self._av_weights = np.zeros(size, float)

        offset = 0
        nsphere = len(self._nlls)

        radii = self._rgrid.radii
        rweights = self._rgrid.weights

        for i in xrange(nsphere):
            nll = self._nlls[i]
            my_points = points[offset:offset+nll]
            my_av_weights = self._av_weights[offset:offset+nll]
            my_weights = weights[offset:offset+nll]

            lebedev_laikov_sphere(my_points, my_av_weights)
            my_points *= radii[i]
            if random_rotate:
                rotmat = get_random_rotation()
                my_points[:] = np.dot(my_points, rotmat)
            my_weights[:] = my_av_weights
            my_weights *= rweights[i]

            offset += nll

        points[:] += center

        IntGrid.__init__(self, points, weights)
        self._log_init()

    def _get_number(self):
        '''The element number of the grid.'''
        return self._number

    number = property(_get_number)

    def _get_center(self):
        '''The center of the grid.'''
        return self._center

    center = property(_get_center)

    def _get_rtransform(self):
        '''The RTransform object of the grid.'''
        return self._rtransform

    rtransform = property(_get_rtransform)

    def _get_int1d(self):
        '''The 1D radial integrator object of the grid.'''
        return self._int1d

    int1d = property(_get_int1d)

    def _get_nlls(self):
        '''The number of Lebedev-Laikov grid points at each sphere.'''
        return self._nlls

    nlls = property(_get_nlls)

    def _get_rgrid(self):
        '''The radial integration grid'''
        return self._rgrid

    rgrid = property(_get_rgrid)

    def _get_nsphere(self):
        '''The number of spheres in the grid.'''
        return len(self._nlls)

    nsphere = property(_get_nsphere)

    def _get_random_rotate(self):
        '''The random rotation flag.'''
        return self._random_rotate

    random_rotate = property(_get_random_rotate)

    def _get_av_weights(self):
        '''The weights needed to compute spherical averages.'''
        return self._av_weights

    av_weights = property(_get_av_weights)

    def _log_init(self):
        if log.do_medium:
            log('Initialized: %s' % self)
            log.deflist([
                ('Size', self.size),
                ('Number of radii', self.nsphere),
                ('Min LL sphere', self._nlls.min()),
                ('Max LL sphere', self._nlls.max()),
                ('Radial Transform', self._rtransform.to_string()),
                ('1D Integrator', self._int1d),
            ])
            # Cite reference
            log.cite('lebedev1999', 'for the use of Lebedev-Laikov grids (quadrature on a sphere)')
            log.blank()

    def get_spherical_average(self, args, output=None):
        '''Returns the spherical average on the radial grid of the product of the given functions'''
        if isinstance(args, np.ndarray) and args.shape == self.shape:
            args = [args]
        args = [arg.ravel() for arg in args if arg is not None]
        args.append(self.av_weights)
        if output is None:
            output = self._rgrid.zeros()
        else:
            assert output.shape == self._rgrid.shape
        # TODO: merge all dot_multi variants in one general-purpose implementation
        dot_multi_parts(args, self._nlls, output)
        return output



# TODO: change this to use the RadialIntGrid object instead of RTransform and Int1D separately.
def interpret_atspec(number, atspec):
    '''Convert atspec to (rtransform, int1d, nlls) tuple

       The atspec argument may be a string refering to a built-in grid file (see
       data/grid) or a tuple with three elements: ``(rtransform, integrator1d,
       nll)`` where:

       * ``rtransform`` is an instance of a subclass of the RTransform class.

       * ``int1d`` is an instance of a subclass of the Integrator1D class.

       * ``nlls`` is a number Lebedev-Laikov grid points for each radial
         grid point. When this argument is not a list, all radial grid
         points get the same Lebedev-Laikov grid

       The number argument is the element number for which the grid is created.
       It is only used when atspec is a string referring to one of the built-in
       grids.
    '''
    if isinstance(atspec, basestring):
        # load
        if atspec not in atgrid_families:
            raise ValueError('Unknown built-in grid: %s' % atspec)
        rtransform, int1d, nlls = atgrid_families[atspec].get(number)
    elif hasattr(atspec, '__iter__') and len(atspec) == 3:
        rtransform, int1d, nlls = atspec

    # Make sure nlls is an array
    if hasattr(nlls, '__iter__'):
        nlls = np.array(nlls, dtype=int)
        if len(nlls) != rtransform.npoint:
            raise ValueError('The size of the radial grid must match the number of elements in nlls')
    else:
        nlls = np.array([nlls]*rtransform.npoint, dtype=int)

    # Finally check the validaty of the nlls
    for nll in nlls:
        if nll not in lebedev_laikov_npoints:
            raise ValueError('A Lebedev-Laikov grid with %i points is not supported.' % nll)

    return rtransform, int1d, nlls


def get_random_rotation():
    '''Return a random rotation matrix'''
    # Get a random unit vector for the axis
    while True:
        axis = np.random.uniform(-1, 1, 3)
        norm = np.linalg.norm(axis)
        if norm < 1.0 and norm > 0.1:
            axis /= norm
            break
    x, y, z = axis

    # Get a random rotation angle
    angle = np.random.uniform(0, 2*np.pi)
    c = np.cos(angle)
    s = np.sin(angle)

    # Rodrigues' rotation formula
    return np.array([
        [x*x*(1-c)+c  , x*y*(1-c)-z*s, x*z*(1-c)+y*s],
        [x*y*(1-c)+z*s, y*y*(1-c)+c  , y*z*(1-c)-x*s],
        [x*z*(1-c)-y*s, y*z*(1-c)+x*s, z*z*(1-c)+c  ],
    ])


class ATGridFamily(object):
    def __init__(self, name):
        self.name = name
        self.members = None

    def get(self, number):
        if self.members is None:
            self._load()

        atspec = self.members.get(number)
        if atspec is None:
            raise ValueError('The atomic grid family %s does not support element %i' % (self.name, number))
        return atspec

    def _load(self):
        fn = context.get_fn('grids/%s.txt' % self.name)
        int1d = SimpsonIntegrator1D()
        self.members = {}
        with open(fn) as f:
            state = 0
            for line in f:
                line = line[:line.find('#')].strip()
                if len(line) > 0:
                    if state == 0:
                        # read element number
                        number = int(line)
                        state = 1
                    elif state == 1:
                        # read rtf string
                        rtf = RTransform.from_string(line)
                        state = 2
                    elif state == 2:
                        nlls = np.array([int(w) for w in line.split()])
                        state = 0
                        self.members[number] = rtf, int1d, nlls


atgrid_families = [
    ATGridFamily('tv-13.1-3'),
    ATGridFamily('tv-13.1-4'),
    ATGridFamily('tv-13.1-5'),
    ATGridFamily('tv-13.1-6'),
]
atgrid_families = dict((af.name.lower(), af) for af in atgrid_families)
