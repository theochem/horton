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
cimport numpy as np
np.import_array()


cdef extern from "xc.h":
    enum: XC_UNPOLARIZED
    enum: XC_POLARIZED

    ctypedef struct xc_func_type:
        pass

    int xc_functional_get_number(char *name)
    bint xc_func_init(xc_func_type *p, int functional, int nspin)
    void xc_lda_end(xc_func_type *p)
    void xc_lda_exc(xc_func_type *p, int npoint, double *rho, double *zk)
    void xc_lda_vxc(xc_func_type *p, int npoint, double *rho, double *vrho)



cdef class LibXCWrapper(object):
    cdef xc_func_type _func_pol
    cdef xc_func_type _func_unpol
    cdef int _func_id

    def __cinit__(self, bytes name):
        '''
           **Arguments:**

           name
                The name of the functional in LibXC
        '''
        self._func_id = -1
        self._func_id = xc_functional_get_number(name)
        if self._func_id < 0:
            raise ValueError('Unknown LibXC functional name: %s' % name)
        retcode = xc_func_init(&self._func_pol, self._func_id, XC_POLARIZED)
        if retcode != 0:
            raise ValueError('Could not initialize polarized LibXC functional: %s' % name)
        retcode = xc_func_init(&self._func_unpol, self._func_id, XC_UNPOLARIZED)
        if retcode != 0:
            raise ValueError('Could not initialize unpolarized LibXC functional: %s' % name)

    def __dealloc__(self):
        if self._func_id >= 0:
            xc_lda_end(&self._func_pol)
            xc_lda_end(&self._func_unpol)

    def compute_exc_unpol(self, np.ndarray[double, ndim=1] rho, np.ndarray[double, ndim=1] zk):
        assert rho.flags['C_CONTIGUOUS']
        npoint = rho.shape[0]
        assert zk.flags['C_CONTIGUOUS']
        assert npoint == zk.shape[0]
        xc_lda_exc(&self._func_unpol, npoint, <double*>rho.data, <double*>zk.data)

    def compute_exc_pol(self, np.ndarray[double, ndim=2] rho, np.ndarray[double, ndim=1] zk):
        assert rho.flags['C_CONTIGUOUS']
        npoint = rho.shape[0]
        assert rho.shape[1] == 2
        assert zk.flags['C_CONTIGUOUS']
        assert npoint == zk.shape[0]
        xc_lda_exc(&self._func_pol, npoint, <double*>rho.data, <double*>zk.data)

    def compute_vxc_unpol(self, np.ndarray[double, ndim=1] rho, np.ndarray[double, ndim=1] vxc):
        assert rho.flags['C_CONTIGUOUS']
        npoint = rho.shape[0]
        assert vxc.flags['C_CONTIGUOUS']
        assert npoint == vxc.shape[0]
        xc_lda_vxc(&self._func_unpol, npoint, <double*>rho.data, <double*>vxc.data)

    def compute_vxc_pol(self, np.ndarray[double, ndim=2] rho, np.ndarray[double, ndim=2] vxc):
        assert rho.flags['C_CONTIGUOUS']
        npoint = rho.shape[0]
        assert rho.shape[1] == 2
        assert vxc.flags['C_CONTIGUOUS']
        assert npoint == vxc.shape[0]
        assert vxc.shape[1] == 2
        xc_lda_vxc(&self._func_pol, npoint, <double*>rho.data, <double*>vxc.data)
