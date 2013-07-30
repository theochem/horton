# -*- coding: utf-8 -*-
# Horton is a development platform for electronic structure methods.
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
#pylint: skip-file


from horton import *
from nose.tools import assert_raises

import numpy as np


def test_consistency():
    for npoint, lmax in lebedev_laikov_npoints.iteritems():
        assert lebedev_laikov_lmaxs[lmax] == npoint


def test_lebedev_laikov_sphere():
    previous_npoint = None
    for i in xrange(1, 132):
        npoint = lebedev_laikov_lmaxs[i]
        if npoint > previous_npoint:
            points = np.zeros((npoint, 3), float)
            weights = np.zeros(npoint, float)
            lebedev_laikov_sphere(points, weights)
            assert abs(weights.sum() - 1.0) < 1e-13
            assert abs(points[:,0].sum()) < 1e-10
            assert abs(points[:,1].sum()) < 1e-10
            assert abs(points[:,2].sum()) < 1e-10
            assert abs(np.dot(points[:,0], weights)) < 1e-15
            assert abs(np.dot(points[:,1], weights)) < 1e-15
            assert abs(np.dot(points[:,2], weights)) < 1e-15
        previous_npoint = npoint
