# -*- coding: utf-8 -*-
# HORTON: Helpful Open-source Research TOol for N-fermion systems.
# Copyright (C) 2011-2019 The HORTON Development Team
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
"""C++ extensions"""


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
    void xc_lda_fxc(xc_func_type *p, int npoint, double *rho, double *v2rho2)
    void xc_gga_exc(xc_func_type *p, int npoint, double *rho, double *sigma, double *zk)
    void xc_gga_vxc(xc_func_type *p, int npoint, double *rho, double *sigma,
                    double *vrho, double *vsigma)
    void xc_gga_fxc(xc_func_type *p, int npoint, double *rho, double *sigma,
                    double *v2rho2, double *v2rhosigma, double *v2sigma2)
    void xc_mgga_exc(xc_func_type *p, int npoint, double *rho, double *sigma,
                     double *lapl, double *tau, double *zk)
    void xc_mgga_vxc(xc_func_type *p, int npoint, double *rho, double *sigma,
                     double *lapl, double *tau, double* vrho, double* vsigma,
                     double* vlapl, double* vtau);
    double xc_hyb_exx_coef(xc_func_type *p)


cdef class LibXCWrapper(object):
    """Base class for restricted and unrestricted LibXC wrappers."""

    cdef xc_func_type _func
    cdef int _func_id
    cdef bytes _key

    def __cinit__(self, bytes key):
        """Initialize a LibXCWrapper.

        Parameters
        ----------
        key : str
            The name of the functional in LibXC, e.g. `"lda_x"`.
        """
        self._key = key
        self._func_id = -1
        self._func_id = xc_functional_get_number(key)
        if self._func_id < 0:
            raise ValueError('Unknown LibXC functional: %s' % key)

    def __dealloc__(self):
        """Deallocate low-level stuff."""
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
        """Return the amount of Hartree-Fock exchange to be used."""
        return xc_hyb_exx_coef(&self._func)



cdef class RLibXCWrapper(LibXCWrapper):
    def __cinit__(self, bytes key):
        """Initialize a RLibXCWrapper.

        Parameters
        ----------
        key : str
            The name of the functional in LibXC, e.g. `"lda_x"`.
        """
        retcode = xc_func_init(&self._func, self._func_id, XC_UNPOLARIZED)
        if retcode != 0:
            self._func_id = -1
            raise ValueError('Could not initialize unpolarized LibXC functional: %s' % key)

    ## LDA

    def compute_lda_exc(self, np.ndarray[double, ndim=1] rho not None,
                              np.ndarray[double, ndim=1] zk not None):
        """Compute the LDA energy density.

        Parameters
        ----------
        rho : np.ndarray, shape=(npoint,)
            The total electron density.
        zk : np.ndarray, shape=(npoint,), output
            The energy density.
        """
        assert rho.flags['C_CONTIGUOUS']
        npoint = rho.shape[0]
        assert zk.flags['C_CONTIGUOUS']
        assert zk.shape[0] == npoint
        xc_lda_exc(&self._func, npoint, &rho[0], &zk[0])

    def compute_lda_vxc(self, np.ndarray[double, ndim=1] rho not None,
                              np.ndarray[double, ndim=1] vrho not None):
        """Compute the LDA potential.

        Parameters
        ----------
        rho : np.ndarray, shape=(npoint,)
            The total electron density.
        vrho : np.ndarray, shape=(npoint,), output
            The LDA potential.
        """
        assert rho.flags['C_CONTIGUOUS']
        npoint = rho.shape[0]
        assert vrho.flags['C_CONTIGUOUS']
        assert vrho.shape[0] == npoint
        xc_lda_vxc(&self._func, npoint, &rho[0], &vrho[0])

    def compute_lda_fxc(self, np.ndarray[double, ndim=1] rho not None,
                              np.ndarray[double, ndim=1] v2rho2 not None):
        """Compute the LDA hardness kernel.

        Parameters
        ----------
        rho : np.ndarray, shape=(npoint,)
            The total electron density.
        v2rho2 : np.ndarray, shape=(npoint,), output
            The (diagonal) LDA kernel.
        """
        assert rho.flags['C_CONTIGUOUS']
        npoint = rho.shape[0]
        assert v2rho2.flags['C_CONTIGUOUS']
        assert v2rho2.shape[0] == npoint
        xc_lda_fxc(&self._func, npoint, &rho[0], &v2rho2[0])

    ## GGA

    def compute_gga_exc(self, np.ndarray[double, ndim=1] rho not None,
                              np.ndarray[double, ndim=1] sigma not None,
                              np.ndarray[double, ndim=1] zk not None):
        """Compute the GGA energy density.

        Parameters
        ----------
        rho : np.ndarray, shape=(npoint,)
            The total electron density.
        sigma : np.ndarray, shape=(npoint,)
            The reduced density gradient norm.
        zk : np.ndarray, shape=(npoint,), output
            The energy density.
        """
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
        """Compute the GGA functional derivatives.

        For every input `x`, a functional derivative is computed, `vx`, stored in an array
        with the same shape.

        Parameters
        ----------
        rho : np.ndarray, shape=(npoint,)
            The total electron density.
        sigma : np.ndarray, shape=(npoint,)
            The reduced density gradient norm.
        vrho : np.ndarray, shape=(npoint,), output
            The LDA part of the potential.
        vsigma : np.ndarray, shape=(npoint,), output
            The GGA part of the potential, i.e. derivatives of the density w.r.t. the
            reduced gradient norm.
        """
        assert rho.flags['C_CONTIGUOUS']
        npoint = rho.shape[0]
        assert sigma.flags['C_CONTIGUOUS']
        assert sigma.shape[0] == npoint
        assert vrho.flags['C_CONTIGUOUS']
        assert vrho.shape[0] == npoint
        assert vsigma.flags['C_CONTIGUOUS']
        assert vsigma.shape[0] == npoint
        xc_gga_vxc(&self._func, npoint, &rho[0], &sigma[0], &vrho[0], &vsigma[0])

    def compute_gga_fxc(self, np.ndarray[double, ndim=1] rho not None,
                              np.ndarray[double, ndim=1] sigma not None,
                              np.ndarray[double, ndim=1] v2rho2 not None,
                              np.ndarray[double, ndim=1] v2rhosigma not None,
                              np.ndarray[double, ndim=1] v2sigma2 not None):
        """Compute the GGA hardness kernel.

        Parameters
        ----------
        rho : np.ndarray, shape=(npoint,)
            The total electron density.
        sigma : np.ndarray, shape=(npoint,)
            The reduced density gradient norm.
        v2rho2 : np.ndarray, shape=(npoint,)
            The second derivative of the energy w.r.t. density (twice).
        v2rhosigma: np.ndarray, shape=(npoint,)
            The second derivative of the energy w.r.t. density (once) and sigma (once).
        v2sigma2: np.ndarray, shape=(npoint,)
            The second derivative of the energy w.r.t. sigma (twice).
        """
        assert rho.flags['C_CONTIGUOUS']
        npoint = rho.shape[0]
        assert sigma.flags['C_CONTIGUOUS']
        assert sigma.shape[0] == npoint
        assert v2rho2.flags['C_CONTIGUOUS']
        assert v2rho2.shape[0] == npoint
        assert v2rhosigma.flags['C_CONTIGUOUS']
        assert v2rhosigma.shape[0] == npoint
        assert v2sigma2.flags['C_CONTIGUOUS']
        assert v2sigma2.shape[0] == npoint
        xc_gga_fxc(&self._func, npoint, &rho[0], &sigma[0], &v2rho2[0], &v2rhosigma[0], &v2sigma2[0])

    ## MGGA

    def compute_mgga_exc(self, np.ndarray[double, ndim=1] rho not None,
                               np.ndarray[double, ndim=1] sigma not None,
                               np.ndarray[double, ndim=1] lapl not None,
                               np.ndarray[double, ndim=1] tau not None,
                               np.ndarray[double, ndim=1] zk not None):
        """Compute the MGGA energy density.

        Parameters
        ----------
        rho : np.ndarray, shape=(npoint,)
            The total electron density.
        sigma : np.ndarray, shape=(npoint,)
            The reduced density gradient norm.
        lapl : np.ndarray, shape=(npoint,)
            The laplacian of the density.
        tau : np.ndarray, shape=(npoint,)
            The kinetic energy density.
        zk : np.ndarray, shape=(npoint,), output
            The energy density.
        """
        assert rho.flags['C_CONTIGUOUS']
        npoint = rho.shape[0]
        assert sigma.flags['C_CONTIGUOUS']
        assert sigma.shape[0] == npoint
        assert lapl.flags['C_CONTIGUOUS']
        assert lapl.shape[0] == npoint
        assert tau.flags['C_CONTIGUOUS']
        assert tau.shape[0] == npoint
        assert zk.flags['C_CONTIGUOUS']
        assert zk.shape[0] == npoint
        xc_mgga_exc(&self._func, npoint, &rho[0], &sigma[0], &lapl[0], &tau[0], &zk[0])

    def compute_mgga_vxc(self, np.ndarray[double, ndim=1] rho not None,
                               np.ndarray[double, ndim=1] sigma not None,
                               np.ndarray[double, ndim=1] lapl not None,
                               np.ndarray[double, ndim=1] tau not None,
                               np.ndarray[double, ndim=1] vrho not None,
                               np.ndarray[double, ndim=1] vsigma not None,
                               np.ndarray[double, ndim=1] vlapl not None,
                               np.ndarray[double, ndim=1] vtau not None,):
        """Compute the MGGA functional derivatives.

        For every input `x`, a functional derivative is computed, `vx`, stored in an array
        with the same shape.

        Parameters
        ----------
        rho : np.ndarray, shape=(npoint,)
            The total electron density.
        sigma : np.ndarray, shape=(npoint,)
            The reduced density gradient norm.
        lapl : np.ndarray, shape=(npoint,)
            The laplacian of the density.
        tau : np.ndarray, shape=(npoint,)
            The kinetic energy density.
        vrho : np.ndarray, shape=(npoint,)
            The derivative of the energy w.r.t. the electron density.
        vsigma : np.ndarray, shape=(npoint,)
            The derivative of the energy w.r.t. the reduced density gradient norm.
        vlapl : np.ndarray, shape=(npoint,)
            The derivative of the energy w.r.t. the laplacian of the density.
        vtau : np.ndarray, shape=(npoint,)
            The derivative of the energy w.r.t. the kinetic energy density.
        """
        assert rho.flags['C_CONTIGUOUS']
        npoint = rho.shape[0]
        assert sigma.flags['C_CONTIGUOUS']
        assert sigma.shape[0] == npoint
        assert lapl.flags['C_CONTIGUOUS']
        assert lapl.shape[0] == npoint
        assert tau.flags['C_CONTIGUOUS']
        assert tau.shape[0] == npoint
        assert vrho.flags['C_CONTIGUOUS']
        assert vrho.shape[0] == npoint
        assert vsigma.flags['C_CONTIGUOUS']
        assert vsigma.shape[0] == npoint
        assert vlapl.flags['C_CONTIGUOUS']
        assert vlapl.shape[0] == npoint
        assert vtau.flags['C_CONTIGUOUS']
        assert vtau.shape[0] == npoint
        xc_mgga_vxc(&self._func, npoint, &rho[0], &sigma[0], &lapl[0], &tau[0], &vrho[0],
                    &vsigma[0], &vlapl[0], &vtau[0])


cdef class ULibXCWrapper(LibXCWrapper):
    def __cinit__(self, bytes key):
        """Initialize a ULibXCWrapper.

        Parameters
        ----------
        key
            The name of the functional in LibXC, e.g. lda_x
        """
        retcode = xc_func_init(&self._func, self._func_id, XC_POLARIZED)
        if retcode != 0:
            self._func_id = -1
            raise ValueError('Could not initialize unpolarized LibXC functional: %s' % key)

    ## LDA

    def compute_lda_exc(self, np.ndarray[double, ndim=2] rho not None,
                              np.ndarray[double, ndim=1] zk not None):
        """Compute the LDA energy density.

        Parameters
        ----------
        rho : np.ndarray, shape=(npoint, 2)
            The alpha and beta electron density.
        zk : np.ndarray, shape=(npoint,), output
            The energy density.
        """
        assert rho.flags['C_CONTIGUOUS']
        npoint = rho.shape[0]
        assert rho.shape[1] == 2
        assert zk.flags['C_CONTIGUOUS']
        assert zk.shape[0] == npoint
        xc_lda_exc(&self._func, npoint, &rho[0, 0], &zk[0])

    def compute_lda_vxc(self, np.ndarray[double, ndim=2] rho not None,
                              np.ndarray[double, ndim=2] vrho not None):
        """Compute the LDA potentials (alpha and beta).

        Parameters
        ----------
        rho : np.ndarray, shape=(npoint, 2)
            The alpha and beta electron density.
        vrho : np.ndarray, shape=(npoint, 2), output
            The SLDA potential, i.e. derivative of the energy w.r.t. alpha and beta
            density.
        """
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
        """Compute the GGA energy density.

        Parameters
        ----------
        rho : np.ndarray, shape=(npoint, 2)
            The alpha and beta electron density.
        sigma : np.ndarray, shape=(npoint, 3)
            The reduced density gradient norms (alpha, alpha), (alpha, beta) and (beta,
            beta).
        zk : np.ndarray, shape=(npoint,), output
            The energy density.
        """
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
        """Compute the GGA functional derivatives.

        For every input `x`, a functional derivative is computed, `vx`, stored in an array
        with the same shape.

        Parameters
        ----------
        rho : np.ndarray, shape=(npoint, 2)
            The alpha and beta electron density.
        sigma : np.ndarray, shape=(npoint, 3)
            The reduced density gradient norms (alpha, alpha), (alpha, beta) and (beta,
            beta).
        vrho : np.ndarray, shape=(npoint, 2), output
            Derivative of the energy w.r.t. alpha and beta density
        vsigma : np.ndarray, shape=(npoint, 3), output
            Derivative of the energy w.r.t reduced density gradient norms.
        """
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
        xc_gga_vxc(&self._func, npoint, &rho[0, 0], &sigma[0, 0], &vrho[0, 0],
                   &vsigma[0, 0])

    ## MGGA

    def compute_mgga_exc(self, np.ndarray[double, ndim=2] rho not None,
                               np.ndarray[double, ndim=2] sigma not None,
                               np.ndarray[double, ndim=2] lapl not None,
                               np.ndarray[double, ndim=2] tau not None,
                               np.ndarray[double, ndim=1] zk not None):
        """Compute the MGGA energy density.

        Parameters
        ----------
        rho : np.ndarray, shape=(npoint, 2)
            The alpha and beta electron density.
        sigma : np.ndarray, shape=(npoint, 3)
            The reduced density gradient norms (alpha, alpha), (alpha, beta) and (beta, beta).
        lapl : np.ndarray, shape=(npoint, 2)
            The laplacian of the alpha and beta density.
        tau : np.ndarray, shape=(npoint, 2)
            The alph and beta kinetic energy density.
        zk : np.ndarray, shape=(npoint,), output
            The energy density.
        """
        assert rho.flags['C_CONTIGUOUS']
        npoint = rho.shape[0]
        assert rho.shape[1] == 2
        assert sigma.flags['C_CONTIGUOUS']
        assert sigma.shape[0] == npoint
        assert sigma.shape[1] == 3
        assert lapl.flags['C_CONTIGUOUS']
        assert lapl.shape[0] == npoint
        assert lapl.shape[1] == 2
        assert tau.flags['C_CONTIGUOUS']
        assert tau.shape[0] == npoint
        assert tau.shape[1] == 2
        assert zk.flags['C_CONTIGUOUS']
        assert zk.shape[0] == npoint
        xc_mgga_exc(&self._func, npoint, &rho[0, 0], &sigma[0, 0], &lapl[0, 0],
                    &tau[0, 0], &zk[0])

    def compute_mgga_vxc(self, np.ndarray[double, ndim=2] rho not None,
                               np.ndarray[double, ndim=2] sigma not None,
                               np.ndarray[double, ndim=2] lapl not None,
                               np.ndarray[double, ndim=2] kin not None,
                               np.ndarray[double, ndim=2] vrho not None,
                               np.ndarray[double, ndim=2] vsigma not None,
                               np.ndarray[double, ndim=2] vlapl not None,
                               np.ndarray[double, ndim=2] vtau not None):
        """Compute the MGGA functional derivatives.

        For every input `x`, a functional derivative is computed, `vx`, stored in an array
        with the same shape.

        Parameters
        ----------
        rho : np.ndarray, shape=(npoint, 2)
            The alpha and beta electron density.
        sigma : np.ndarray, shape=(npoint, 3)
            The reduced density gradient norms (alpha, alpha), (alpha, beta) and (beta, beta).
        lapl : np.ndarray, shape=(npoint, 2)
            The laplacian of the alpha and beta density.
        tau : np.ndarray, shape=(npoint, 2)
            The alph and beta kinetic energy density.
        vrho : np.ndarray, shape=(npoint, 2)
            The derivative of the energy w.r.t. the alpha and beta electron density.
        vsigma : np.ndarray, shape=(npoint, 3)
            The derivative of the energy w.r.t. the reduced density gradient norms.
        vlapl : np.ndarray, shape=(npoint, 2)
            The derivative of the energy w.r.t. the laplacian of the alpha and beta density.
        vtau : np.ndarray, shape=(npoint, 2)
            The derivative of the energy w.r.t. the alpha and beta kinetic energy density.
        """
        assert rho.flags['C_CONTIGUOUS']
        npoint = rho.shape[0]
        assert rho.shape[1] == 2
        assert sigma.flags['C_CONTIGUOUS']
        assert sigma.shape[0] == npoint
        assert sigma.shape[1] == 3
        assert lapl.flags['C_CONTIGUOUS']
        assert lapl.shape[0] == npoint
        assert lapl.shape[1] == 2
        assert kin.flags['C_CONTIGUOUS']
        assert kin.shape[0] == npoint
        assert kin.shape[1] == 2
        assert vrho.flags['C_CONTIGUOUS']
        assert vrho.shape[0] == npoint
        assert vrho.shape[1] == 2
        assert vsigma.flags['C_CONTIGUOUS']
        assert vsigma.shape[0] == npoint
        assert vsigma.shape[1] == 3
        assert vlapl.flags['C_CONTIGUOUS']
        assert vlapl.shape[0] == npoint
        assert vlapl.shape[1] == 2
        assert vtau.flags['C_CONTIGUOUS']
        assert vtau.shape[0] == npoint
        assert vtau.shape[1] == 2
        xc_mgga_vxc(&self._func, npoint, &rho[0, 0], &sigma[0, 0], &lapl[0, 0],
                    &kin[0, 0], &vrho[0, 0], &vsigma[0, 0], &vlapl[0, 0], &vtau[0, 0])
