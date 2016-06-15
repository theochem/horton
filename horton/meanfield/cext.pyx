# -*- coding: utf-8 -*-
# HORTON: Helpful Open-source Research TOol for N-fermion systems.
# Copyright (C) 2011-2016 The HORTON Development Team
#
# This file is part of HORTON.
#
# HORTON is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
#
# HORTON is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, see <http://www.gnu.org/licenses/>
#
# --
'''C++ extensions'''


import numpy as np
cimport numpy as np
np.import_array()


__all__ = [
    'RLibXCWrapper', 'ULibXCWrapper'
]


cdef extern from "xc.h":
    enum: XC_UNPOLARIZED
    enum: XC_POLARIZED

    ctypedef struct xc_func_info_type:
        int number
        int kind
        char* name
        int family
        char* refs

    ctypedef struct xc_func_type:
        xc_func_info_type* info

    int xc_functional_get_number(char *name)
    bint xc_func_init(xc_func_type *p, int functional, int nspin)
    void xc_func_end(xc_func_type *p)
    void xc_lda_exc(xc_func_type *p, int npoint, double *rho, double *zk)
    void xc_lda_vxc(xc_func_type *p, int npoint, double *rho, double *vrho)
    void xc_gga_exc(xc_func_type *p, int npoint, double *rho, double *sigma, double *zk)
    void xc_gga_vxc(xc_func_type *p, int npoint, double *rho, double *sigma, double *vrho, double *vsigma)
    double xc_hyb_exx_coef(xc_func_type *p)


cdef class LibXCWrapper(object):
    cdef xc_func_type _func
    cdef int _func_id
    cdef bytes _key

    def __cinit__(self, bytes key):
        '''
           **Arguments:**

           key
                The name of the functional in LibXC, e.g. lda_x
        '''
        self._key = key
        self._func_id = -1
        self._func_id = xc_functional_get_number(key)
        if self._func_id < 0:
            raise ValueError('Unknown LibXC functional: %s' % key)

    def __dealloc__(self):
        if self._func_id >= 0:
            xc_func_end(&self._func)

    ## INFO

    property key:
        def __get__(self):
            return self._key

    property number:
        def __get__(self):
            return self._func_id

    property kind:
        def __get__(self):
            return self._func.info[0].kind

    property name:
        def __get__(self):
            return self._func.info[0].name

    property family:
        def __get__(self):
            return self._func.info[0].family

    property refs:
        def __get__(self):
            return self._func.info[0].refs

    ## HYB GGA

    def get_hyb_exx_fraction(self):
        return xc_hyb_exx_coef(&self._func)



cdef class RLibXCWrapper(LibXCWrapper):
    def __cinit__(self, bytes key):
        '''
           **Arguments:**

           key
                The name of the functional in LibXC, e.g. lda_x
        '''
        retcode = xc_func_init(&self._func, self._func_id, XC_UNPOLARIZED)
        if retcode != 0:
            self._func_id = -1
            raise ValueError('Could not initialize unpolarized LibXC functional: %s' % key)

    ## LDA

    def compute_lda_exc(self, np.ndarray[double, ndim=1] rho not None,
                              np.ndarray[double, ndim=1] zk not None):
        assert rho.flags['C_CONTIGUOUS']
        npoint = rho.shape[0]
        assert zk.flags['C_CONTIGUOUS']
        assert zk.shape[0] == npoint
        xc_lda_exc(&self._func, npoint, &rho[0], &zk[0])

    def compute_lda_vxc(self, np.ndarray[double, ndim=1] rho not None,
                              np.ndarray[double, ndim=1] vrho not None):
        assert rho.flags['C_CONTIGUOUS']
        npoint = rho.shape[0]
        assert vrho.flags['C_CONTIGUOUS']
        assert vrho.shape[0] == npoint
        xc_lda_vxc(&self._func, npoint, &rho[0], &vrho[0])

    ## GGA

    def compute_gga_exc(self, np.ndarray[double, ndim=1] rho not None,
                              np.ndarray[double, ndim=1] sigma not None,
                              np.ndarray[double, ndim=1] zk not None):
        assert rho.flags['C_CONTIGUOUS']
        npoint = rho.shape[0]
        assert sigma.flags['C_CONTIGUOUS']
        assert sigma.shape[0] == npoint
        assert zk.flags['C_CONTIGUOUS']
        assert zk.shape[0] == npoint
        xc_gga_exc(&self._func, npoint, &rho[0], &sigma[0], &zk[0])

    def compute_gga_vxc(self, np.ndarray[double, ndim=1] rho not None,
                              np.ndarray[double, ndim=1] sigma not None,
                              np.ndarray[double, ndim=1] vrho not None,
                              np.ndarray[double, ndim=1] vsigma not None):
        assert rho.flags['C_CONTIGUOUS']
        npoint = rho.shape[0]
        assert sigma.flags['C_CONTIGUOUS']
        assert sigma.shape[0] == npoint
        assert vrho.flags['C_CONTIGUOUS']
        assert vrho.shape[0] == npoint
        assert vsigma.flags['C_CONTIGUOUS']
        assert vsigma.shape[0] == npoint
        xc_gga_vxc(&self._func, npoint, &rho[0], &sigma[0], &vrho[0], &vsigma[0])


cdef class ULibXCWrapper(LibXCWrapper):
    def __cinit__(self, bytes key):
        '''
           **Arguments:**

           key
                The name of the functional in LibXC, e.g. lda_x
        '''
        retcode = xc_func_init(&self._func, self._func_id, XC_POLARIZED)
        if retcode != 0:
            self._func_id = -1
            raise ValueError('Could not initialize unpolarized LibXC functional: %s' % key)

    ## LDA

    def compute_lda_exc(self, np.ndarray[double, ndim=2] rho not None,
                              np.ndarray[double, ndim=1] zk not None):
        assert rho.flags['C_CONTIGUOUS']
        npoint = rho.shape[0]
        assert rho.shape[1] == 2
        assert zk.flags['C_CONTIGUOUS']
        assert zk.shape[0] == npoint
        xc_lda_exc(&self._func, npoint, &rho[0, 0], &zk[0])

    def compute_lda_vxc(self, np.ndarray[double, ndim=2] rho not None,
                              np.ndarray[double, ndim=2] vrho not None):
        assert rho.flags['C_CONTIGUOUS']
        npoint = rho.shape[0]
        assert rho.shape[1] == 2
        assert vrho.flags['C_CONTIGUOUS']
        assert vrho.shape[0] == npoint
        assert vrho.shape[1] == 2
        xc_lda_vxc(&self._func, npoint, &rho[0, 0], &vrho[0, 0])

    ## GGA

    def compute_gga_exc(self, np.ndarray[double, ndim=2] rho not None,
                              np.ndarray[double, ndim=2] sigma not None,
                              np.ndarray[double, ndim=1] zk not None):
        assert rho.flags['C_CONTIGUOUS']
        npoint = rho.shape[0]
        assert rho.shape[1] == 2
        assert sigma.flags['C_CONTIGUOUS']
        assert sigma.shape[1] == 3
        assert sigma.shape[0] == npoint
        assert zk.flags['C_CONTIGUOUS']
        assert zk.shape[0] == npoint
        xc_gga_exc(&self._func, npoint, &rho[0, 0], &sigma[0, 0], &zk[0])

    def compute_gga_vxc(self, np.ndarray[double, ndim=2] rho not None,
                              np.ndarray[double, ndim=2] sigma not None,
                              np.ndarray[double, ndim=2] vrho not None,
                              np.ndarray[double, ndim=2] vsigma not None):
        assert rho.flags['C_CONTIGUOUS']
        npoint = rho.shape[0]
        assert rho.shape[1] == 2
        assert sigma.flags['C_CONTIGUOUS']
        assert sigma.shape[0] == npoint
        assert sigma.shape[1] == 3
        assert vrho.flags['C_CONTIGUOUS']
        assert vrho.shape[0] == npoint
        assert vrho.shape[1] == 2
        assert vsigma.flags['C_CONTIGUOUS']
        assert vsigma.shape[0] == npoint
        assert vsigma.shape[1] == 3
        xc_gga_vxc(&self._func, npoint, &rho[0, 0], &sigma[0, 0], &vrho[0, 0], &vsigma[0, 0])
