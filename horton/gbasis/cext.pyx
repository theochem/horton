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

cimport boys
cimport cartpure
cimport common
cimport gbasis
cimport ints
cimport iter_gb
cimport iter_pow

import atexit


__all__ = [
    # boys
    'boys_function',
    # cartpure
    'cart_to_pure_low',
    # common
    'fac', 'fac2', 'binom', 'get_shell_nbasis', 'get_max_shell_type',
    # gbasis
    'GOBasis',
    # ints
    'gpt_coeff', 'gb_overlap_int1d', 'GB2OverlapIntegral', 'GB2KineticIntegral',
    'nuclear_attraction_helper', 'GB2NuclearAttractionIntegral',
    'GB4ElectronReuplsionIntegralLibInt',
    'gob_normalization',
    # iter_gb
    'IterGB2', 'IterGB4',
    # iter_pow
    'iter_pow1_inc', 'IterPow2',
]


#
# boys wrappers (for testing only)
#


def boys_function(long m, double t):
    return boys.boys_function(m, t)


#
# cartpure wrappers (for testing only)
#


def cart_to_pure_low(np.ndarray[double] work_cart,
                     np.ndarray[double] work_pure, long shell_type,
                     long nant, long npost):
    assert work_cart.flags['C_CONTIGUOUS']
    assert work_pure.flags['C_CONTIGUOUS']
    cartpure.cart_to_pure_low(
        <double*> work_cart.data, <double*> work_pure.data, shell_type, nant,
        npost
    )


#
# common wrappers
#


def fac(long n):
    return common.fac(n)


def fac2(long n):
    return common.fac2(n)


def binom(long n, long m):
    return common.binom(n, m)


def get_shell_nbasis(long shell_type):
    result = common.get_shell_nbasis(shell_type)
    if result <= 0:
        raise ValueError("shell_type -1 is not supported.")
    return result


def get_max_shell_type():
    return common.get_max_shell_type()

#
# gbasis wrappers
#


cdef class GBasis:
    """
       This class describes basis sets applied to a certain molecular structure.

       **Arguments:**

       centers
            A numpy array with centers for the basis functions.
            shape = (ncenter, 3)

       shell_map
            An array with the center index for each shell.
            shape = (nshell,)

       nprims
            The number of primitives in each shell.
            shape = (nshell,)

       shell_types
            An array with contraction types: 0 = S, 1 = P, 2 = Cartesian D,
            3 = Cartesian F, ..., -2 = pure D, -3 = pure F, ...
            shape = (nshell,)

       alphas
            The exponents of the primitives in one shell.
            shape = (sum(nprims),)

       con_coeffs
            The contraction coefficients of the primitives for each
            contraction in a contiguous array. The coefficients are ordered
            according to the shells. Within each shell, the coefficients are
            grouped per exponent.
            shape = (sum(nprims),)

       All arguments may also be lists. In any case, copies are made of the
       arguments and stored internally, and are not meant to be modified once
       the GOBasis object is created.

       The order of the pure shells is based on the order of real spherical.
       The functions are sorted from low to high magnetic quantum number,
       with cosine-like functions before the sine-like functions. The order
       of functions in a Cartesian shell is alhpabetic. Some examples:

       shell_type = 0, S:
         0 -> 1
       shell_type = 1, P:
         0 -> x
         1 -> y
         2 -> z
       shell_type = 2, Cartesian D:
         0 -> xx
         1 -> xy
         2 -> xz
         3 -> yy
         4 -> yz
         5 -> zz
       shell_type = 3, Cartesian F:
         0 -> xxx
         1 -> xxy
         2 -> xxz
         3 -> xyy
         4 -> xyz
         5 -> xzz
         6 -> yyy
         7 -> yyz
         8 -> yzz
         9 -> zzz
       shell_type = -1, not allowed
       shell_type = -2, pure D:
         0 -> zz
         1 -> xz
         2 -> yz
         3 -> xx-yy
         4 -> xy
       shell_type = -3, pure F:
         0 -> zzz
         1 -> xzz
         2 -> yzz
         3 -> xxz-yyz
         4 -> xyz
         5 -> xxx-3xyy
         6 -> 3xxy-yyy
    """

    cdef gbasis.GBasis* _this
    # Keep reference to arrays to make sure they will not be deallocated.
    cdef np.ndarray _centers
    cdef np.ndarray _shell_map
    cdef np.ndarray _nprims
    cdef np.ndarray _shell_types
    cdef np.ndarray _alphas
    cdef np.ndarray _con_coeffs

    def __cinit__(self, centers, shell_map, nprims, shell_types, alphas, con_coeffs):
        # Make private copies of the input arrays.
        self._centers = np.array(centers, dtype=float)
        self._shell_map = np.array(shell_map, dtype=int)
        self._nprims = np.array(nprims, dtype=int)
        self._shell_types = np.array(shell_types, dtype=int)
        self._alphas = np.array(alphas, dtype=float)
        self._con_coeffs = np.array(con_coeffs, dtype=float)

        # Set arrays unwritable because:
        #   (i) derived properties will be stored
        #   (ii) consistency tests below are only performed once (below)
        # In principle, one can make the arrays writable again, but then one
        # is deliberately taking risks. A strict read-only buffer would be
        # ideal but is not possible. One can always pass a pointer to a
        # C library and start messing things up.
        self._centers.flags.writeable = False
        self._shell_map.flags.writeable = False
        self._nprims.flags.writeable = False
        self._shell_types.flags.writeable = False
        self._alphas.flags.writeable = False
        self._con_coeffs.flags.writeable = False

        # Check array dimensions
        if self._centers.ndim != 2:
            raise TypeError('centers must be a 2D array')
        if self._shell_map.ndim != 1:
            raise TypeError('shell_map must be a 1D array')
        if self._nprims.ndim != 1:
            raise TypeError('nprims must be a 1D array')
        if self._shell_types.ndim != 1:
            raise TypeError('shell_types must be a 1D array')
        if self._alphas.ndim != 1:
            raise TypeError('alphas must be a 1D array')
        if self._con_coeffs.ndim != 1:
            raise TypeError('con_coeffs must be a 1D array')

        # Essential array checks
        if self._centers.shape[1] != 3:
            raise TypeError('centers must have three columns.')
        if self._nprims.shape[0] != self._shell_map.shape[0]:
            raise TypeError('nprims and shell_map must have the same length.')
        if self._shell_types.shape[0] != self._shell_map.shape[0]:
            raise TypeError('shell_types and shell_map must have the same length.')
        if self._alphas.shape[0] != self._con_coeffs.shape[0]:
            raise TypeError('alphas and con_coeffs must have the same length.')

        # Consistency checks
        if self._shell_map.min() < 0:
            raise ValueError('shell_map can not contain negative values.')
        if self._shell_map.max() >= self.centers.shape[0]:
            raise ValueError('shell_map can not contain values larger than the number of centers.')
        if self._nprims.min() < 1:
            raise ValueError('nprims elements must be strictly positive.')
        if (self._shell_types == -1).any():
            raise ValueError('The shell_type -1 is not supported.')
        cdef long nprim_total = self._nprims.sum()
        if (self._alphas.shape[0] != nprim_total):
            raise TypeError('The length of alphas must equal the total number of primitives.')

    def to_hdf5(self, grp):
        raise NotImplementedError

    # Array properties

    property centers:
        def __get__(self):
            return self._centers.view()

    property shell_map:
        def __get__(self):
            return self._shell_map.view()

    property nprims:
        def __get__(self):
            return self._nprims.view()

    property shell_types:
        def __get__(self):
            return self._shell_types.view()

    property alphas:
        def __get__(self):
            return self._alphas.view()

    property con_coeffs:
        def __get__(self):
            return self._con_coeffs.view()

    # Array sizes

    property ncenters:
        def __get__(self):
            return self.center.shape[0]

    property nshell:
        def __get__(self):
            return self.shell_map.shape[0]

    property nprim_total:
        def __get__(self):
            return self.nprim_total.shape[0]

    # Other properties

    property nbasis:
        def __get__(self):
            return self._this.get_nbasis()

    property nscales:
        def __get__(self):
            return self._this.get_nscales()

    property max_shell_type:
        def __get__(self):
            return self._this.get_max_shell_type()

    def get_scales(self):
        cdef np.npy_intp shape[1]
        shape[0] = self.nscales
        tmp = np.PyArray_SimpleNewFromData(1, shape, np.NPY_DOUBLE, <void*> self._this.get_scales(0))
        return tmp.copy()


cdef class GOBasis(GBasis):
    def __cinit__(self, centers, shell_map, nprims, shell_types, alphas, con_coeffs):
        self._this = <gbasis.GBasis*> new gbasis.GOBasis(
            <double*>self._centers.data, <long*>self._shell_map.data,
            <long*>self._nprims.data, <long*>self._shell_types.data,
            <double*>self._alphas.data, <double*>self._con_coeffs.data,
            self._centers.shape[0], self._shell_types.shape[0], self._alphas.shape[0]
        )

    @classmethod
    def from_hdf5(cls, grp, lf):
        return GOBasis(
            np.array(grp['centers']),
            np.array(grp['shell_map']),
            np.array(grp['nprims']),
            np.array(grp['shell_types']),
            np.array(grp['alphas']),
            np.array(grp['con_coeffs'])
        )

    def to_hdf5(self, grp):
        grp['centers'] = self.centers
        grp['shell_map'] = self.shell_map
        grp['nprims'] = self.nprims
        grp['shell_types'] = self.shell_types
        grp['alphas'] = self.alphas
        grp['con_coeffs'] = self.con_coeffs

    # Integral computation

    def check_matrix_one_body(self, matrix):
        assert matrix.ndim == 2
        assert matrix.flags['C_CONTIGUOUS']
        assert matrix.shape[0] == self.nbasis
        assert matrix.shape[1] == self.nbasis

    def check_matrix_two_body(self, matrix):
        assert matrix.ndim == 4
        assert matrix.flags['C_CONTIGUOUS']
        assert matrix.shape[0] == self.nbasis
        assert matrix.shape[1] == self.nbasis
        assert matrix.shape[2] == self.nbasis
        assert matrix.shape[3] == self.nbasis

    def compute_overlap(self, overlap):
        """Compute the overlap matrix in a Gaussian orbital basis."""
        cdef np.ndarray output = overlap._array
        self.check_matrix_one_body(output)
        gobasis = <gbasis.GOBasis*>self._this
        (<gbasis.GOBasis*>self._this).compute_overlap(<double*>output.data)

    def compute_kinetic(self, kinetic):
        """Compute the kinetic energy matrix in a Gaussian orbital basis."""
        cdef np.ndarray output = kinetic._array
        self.check_matrix_one_body(output)
        (<gbasis.GOBasis*>self._this).compute_kinetic(<double*>output.data)

    def compute_nuclear_attraction(self, np.ndarray[double, ndim=1] charges,
                                   np.ndarray[double, ndim=2] centers,
                                   nuclear_attraction):
        """Compute the kintic energy matrix in a Gaussian orbital basis."""
        cdef np.ndarray output = nuclear_attraction._array
        self.check_matrix_one_body(output)
        assert charges.flags['C_CONTIGUOUS']
        cdef long ncharge = charges.shape[0]
        assert centers.flags['C_CONTIGUOUS']
        assert centers.shape[0] == ncharge
        (<gbasis.GOBasis*>self._this).compute_nuclear_attraction(
            <double*>charges.data, <double*>centers.data, ncharge,
            <double*>output.data
        )

    def compute_electron_repulsion(self, electron_repulsion):
        cdef np.ndarray output = electron_repulsion._array
        self.check_matrix_two_body(output)
        (<gbasis.GOBasis*>self._this).compute_electron_repulsion(<double*>output.data)


#
# ints wrappers (for testing only)
#


def gpt_coeff(long k, long n0, long n1, double pa, double pb):
    return ints.gpt_coeff(k, n0, n1, pa, pb)


def gb_overlap_int1d(long n0, long n1, double pa, double pb, double inv_gamma):
    return ints.gb_overlap_int1d(n0, n1, pa, pb, inv_gamma)


def gob_normalization(double alpha, np.ndarray[long, ndim=1] n):
    assert n.flags['C_CONTIGUOUS']
    assert n.shape[0] == 3
    return ints.gob_normalization(alpha, <long*>n.data)


cdef class GB2Integral:
    '''Wrapper for ints.GB2Integral, for testing only'''
    cdef ints.GB2Integral* _this

    def __dealloc__(self):
        del self._this

    property nwork:
        def __get__(self):
            return self._this.get_nwork()

    property max_shell_type:
        def __get__(self):
            return self._this.get_max_shell_type()

    property max_nbasis:
        def __get__(self):
            return self._this.get_max_nbasis()

    def reset(self, long shell_type0, long shell_type1,
              np.ndarray[double, ndim=1] r0, np.ndarray[double, ndim=1] r1):
        assert r0.flags['C_CONTIGUOUS']
        assert r0.shape[0] == 3
        assert r1.flags['C_CONTIGUOUS']
        assert r1.shape[0] == 3
        self._this.reset(shell_type0, shell_type1, <double*>r0.data, <double*>r1.data)

    def add(self, double coeff, double alpha0, double alpha1,
            np.ndarray[double, ndim=1] scales0, np.ndarray[double, ndim=1] scales1):
        assert scales0.flags['C_CONTIGUOUS']
        assert scales0.shape[0] == get_shell_nbasis(abs(self._this.get_shell_type0()))
        assert scales1.flags['C_CONTIGUOUS']
        assert scales1.shape[0] == get_shell_nbasis(abs(self._this.get_shell_type1()))
        self._this.add(coeff, alpha0, alpha1, <double*>scales0.data, <double*>scales1.data)

    def cart_to_pure(self):
        self._this.cart_to_pure()

    def get_work(self, shape0, shape1):
        '''This returns a **copy** of the c++ work array.

           Returning a numpy array with a buffer created in c++ is dangerous.
           If the c++ array becomes deallocated, the numpy array may still
           point to the deallocated memory. For that reason, a copy is returned.
           Speed is not an issue as this class is only used for testing.
        '''
        cdef np.npy_intp shape[2]
        assert shape0 > 0
        assert shape1 > 0
        assert shape0 <= self.max_nbasis
        assert shape1 <= self.max_nbasis
        shape[0] = shape0
        shape[1] = shape1
        tmp = np.PyArray_SimpleNewFromData(2, shape, np.NPY_DOUBLE, <void*> self._this.get_work())
        return tmp.copy()


cdef class GB2OverlapIntegral(GB2Integral):
    '''Wrapper for ints.GB2OverlapIntegral, for testing only'''

    def __cinit__(self, long max_nbasis):
        self._this = <ints.GB2Integral*>(new ints.GB2OverlapIntegral(max_nbasis))


cdef class GB2KineticIntegral(GB2Integral):
    '''Wrapper for ints.GB2OverlapIntegral, for testing only'''

    def __cinit__(self, long max_nbasis):
        self._this = <ints.GB2Integral*>(new ints.GB2KineticIntegral(max_nbasis))


def nuclear_attraction_helper(np.ndarray[double, ndim=1] work_g, long n0,
                              long n1, double pa, double pb, double cp,
                              double gamma_inv):
    assert work_g.flags['C_CONTIGUOUS']
    assert work_g.shape[0] == n0+n1+1
    ints.nuclear_attraction_helper(<double*>work_g.data, n0, n1, pa, pb, cp, gamma_inv)


cdef class GB2NuclearAttractionIntegral(GB2Integral):
    '''Wrapper for ints.GB2NuclearAttractionIntegral, for testing only'''
    # make an additional reference to these arguments to avoid deallocation
    cdef np.ndarray _charges
    cdef np.ndarray _centers

    def __cinit__(self, long max_nbasis, np.ndarray[double, ndim=1] charges, np.ndarray[double, ndim=2] centers):
        assert charges.flags['C_CONTIGUOUS']
        cdef long ncharge = charges.shape[0]
        assert centers.flags['C_CONTIGUOUS']
        assert centers.shape[0] == ncharge
        self._charges = charges
        self._centers = centers
        self._this = <ints.GB2Integral*>(new ints.GB2NuclearAttractionIntegral(
            max_nbasis, <double*>charges.data, <double*>centers.data, ncharge
        ))


ints.libint2_static_init()
def libint2_static_cleanup():
    ints.libint2_static_cleanup()
atexit.register(libint2_static_cleanup)


cdef class GB4Integral:
    '''Wrapper for ints.GB4Integral, for testing only'''
    cdef ints.GB4Integral* _this

    def __dealloc__(self):
        del self._this

    property nwork:
        def __get__(self):
            return self._this.get_nwork()

    property max_shell_type:
        def __get__(self):
            return self._this.get_max_shell_type()

    property max_nbasis:
        def __get__(self):
            return self._this.get_max_nbasis()

    def reset(self, long shell_type0, long shell_type1, long shell_type2, long shell_type3,
              np.ndarray[double, ndim=1] r0, np.ndarray[double, ndim=1] r1,
              np.ndarray[double, ndim=1] r2, np.ndarray[double, ndim=1] r3):
        assert r0.flags['C_CONTIGUOUS']
        assert r0.shape[0] == 3
        assert r1.flags['C_CONTIGUOUS']
        assert r1.shape[0] == 3
        assert r2.flags['C_CONTIGUOUS']
        assert r2.shape[0] == 3
        assert r3.flags['C_CONTIGUOUS']
        assert r3.shape[0] == 3
        self._this.reset(shell_type0, shell_type1, shell_type2, shell_type3,
                         <double*>r0.data, <double*>r1.data, <double*>r2.data, <double*>r3.data)

    def add(self, double coeff, double alpha0, double alpha1, double alpha2, double alpha3,
            np.ndarray[double, ndim=1] scales0, np.ndarray[double, ndim=1] scales1,
            np.ndarray[double, ndim=1] scales2, np.ndarray[double, ndim=1] scales3):
        assert scales0.flags['C_CONTIGUOUS']
        assert scales0.shape[0] == get_shell_nbasis(abs(self._this.get_shell_type0()))
        assert scales1.flags['C_CONTIGUOUS']
        assert scales1.shape[0] == get_shell_nbasis(abs(self._this.get_shell_type1()))
        assert scales2.flags['C_CONTIGUOUS']
        assert scales2.shape[0] == get_shell_nbasis(abs(self._this.get_shell_type2()))
        assert scales3.flags['C_CONTIGUOUS']
        assert scales3.shape[0] == get_shell_nbasis(abs(self._this.get_shell_type3()))
        self._this.add(coeff, alpha0, alpha1, alpha2, alpha3,
                       <double*>scales0.data, <double*>scales1.data,
                       <double*>scales2.data, <double*>scales3.data)

    def cart_to_pure(self):
        self._this.cart_to_pure()

    def get_work(self, shape0, shape1, shape2, shape3):
        '''This returns a **copy** of the c++ work array.

           Returning a numpy array with a buffer created in c++ is dangerous.
           If the c++ array becomes deallocated, the numpy array may still
           point to the deallocated memory. For that reason, a copy is returned.
           Speed is not an issue as this class is only used for testing.
        '''
        cdef np.npy_intp shape[4]
        assert shape0 > 0
        assert shape1 > 0
        assert shape2 > 0
        assert shape3 > 0
        assert shape0 <= self.max_nbasis
        assert shape1 <= self.max_nbasis
        assert shape2 <= self.max_nbasis
        assert shape3 <= self.max_nbasis
        shape[0] = shape0
        shape[1] = shape1
        shape[2] = shape2
        shape[3] = shape3
        tmp = np.PyArray_SimpleNewFromData(4, shape, np.NPY_DOUBLE, <void*> self._this.get_work())
        return tmp.copy()


cdef class GB4ElectronReuplsionIntegralLibInt(GB4Integral):
    '''Wrapper for ints.GB4ElectronReuplsionIntegralLibInt, for testing only'''

    def __cinit__(self, long max_nbasis):
        self._this = <ints.GB4Integral*>(new ints.GB4ElectronReuplsionIntegralLibInt(max_nbasis))



#
# iter_gb wrappers (for testing only)
#


cdef class IterGB2:
    """Wrapper for the IterGB2 class, for testing only."""
    cdef iter_gb.IterGB2* _this
    cdef GBasis _gbasis

    def __cinit__(self, GBasis gbasis):
        self._this = new iter_gb.IterGB2(gbasis._this)
        self._gbasis = gbasis

    def __dealloc__(self):
        del self._this

    def inc_shell(self):
        return self._this.inc_shell()

    def update_shell(self):
        self._this.update_shell()

    def inc_prim(self):
        return self._this.inc_prim()

    def update_prim(self):
        self._this.update_prim()

    def store(self, np.ndarray[double, ndim=2] work,
              np.ndarray[double, ndim=2] output):
        max_shell_nbasis = get_shell_nbasis(self._gbasis.max_shell_type)
        assert work.shape[0] == get_shell_nbasis(self._this.shell_type0)
        assert work.shape[1] == get_shell_nbasis(self._this.shell_type1)
        assert work.flags['C_CONTIGUOUS']
        assert output.shape[0] == self._gbasis.nbasis
        assert output.shape[1] == self._gbasis.nbasis
        assert output.flags['C_CONTIGUOUS']
        self._this.store(<double*>work.data, <double*>output.data)

    property public_fields:
        def __get__(self):
            return (
                self._this.con_coeff,
                self._this.shell_type0, self._this.shell_type1,
                self._this.alpha0, self._this.alpha1,
                self._this.r0[0], self._this.r0[1], self._this.r0[2],
                self._this.r1[0], self._this.r1[1], self._this.r1[2],
                self._this.ibasis0, self._this.ibasis1,
            )

    property private_fields:
        def __get__(self):
            return (
                self._this.ishell0, self._this.ishell1,
                self._this.nprim0, self._this.nprim1,
                self._this.oprim0, self._this.oprim1,
                self._this.iprim0, self._this.iprim1,
            )


cdef class IterGB4:
    """Wrapper for the IterGB4 class, for testing only."""
    cdef iter_gb.IterGB4* _this
    cdef GBasis _gbasis

    def __cinit__(self, GBasis gbasis):
        self._this = new iter_gb.IterGB4(gbasis._this)
        self._gbasis = gbasis

    def __dealloc__(self):
        del self._this

    def inc_shell(self):
        return self._this.inc_shell()

    def update_shell(self):
        self._this.update_shell()

    def inc_prim(self):
        return self._this.inc_prim()

    def update_prim(self):
        self._this.update_prim()

    def store(self, np.ndarray[double, ndim=4] work,
              np.ndarray[double, ndim=4] output):
        max_shell_nbasis = get_shell_nbasis(self._gbasis.max_shell_type)
        assert work.shape[0] == get_shell_nbasis(self._this.shell_type0)
        assert work.shape[1] == get_shell_nbasis(self._this.shell_type1)
        assert work.shape[2] == get_shell_nbasis(self._this.shell_type2)
        assert work.shape[3] == get_shell_nbasis(self._this.shell_type3)
        assert work.flags['C_CONTIGUOUS']
        assert output.shape[0] == self._gbasis.nbasis
        assert output.shape[1] == self._gbasis.nbasis
        assert output.shape[2] == self._gbasis.nbasis
        assert output.shape[3] == self._gbasis.nbasis
        assert output.flags['C_CONTIGUOUS']
        self._this.store(<double*>work.data, <double*>output.data)

    property public_fields:
        def __get__(self):
            return (
                self._this.con_coeff,
                self._this.shell_type0, self._this.shell_type1, self._this.shell_type2, self._this.shell_type3,
                self._this.alpha0, self._this.alpha1, self._this.alpha2, self._this.alpha3,
                self._this.r0[0], self._this.r0[1], self._this.r0[2],
                self._this.r1[0], self._this.r1[1], self._this.r1[2],
                self._this.r2[0], self._this.r2[1], self._this.r2[2],
                self._this.r3[0], self._this.r3[1], self._this.r3[2],
                self._this.ibasis0, self._this.ibasis1, self._this.ibasis2, self._this.ibasis3,
            )

    property private_fields:
        def __get__(self):
            return (
                self._this.ishell0, self._this.ishell1, self._this.ishell2, self._this.ishell3,
                self._this.nprim0, self._this.nprim1, self._this.nprim2, self._this.nprim3,
                self._this.oprim0, self._this.oprim1, self._this.oprim2, self._this.oprim3,
                self._this.iprim0, self._this.iprim1, self._this.iprim2, self._this.iprim3,
            )


#
# iter_pow wrappers (for testing only)
#


def iter_pow1_inc(np.ndarray[long, ndim=1] n):
    assert n.flags['C_CONTIGUOUS']
    assert n.shape[0] == 3
    return iter_pow.iter_pow1_inc(<long*>n.data)


cdef class IterPow2:
    """Wrapper for the IterPow2 class, for testing only."""
    cdef iter_pow.IterPow2* _c_i2p

    def __cinit__(self):
        self._c_i2p = new iter_pow.IterPow2()

    def __dealloc__(self):
        del self._c_i2p

    def __init__(self, long shell_type0, long shell_type1, max_nbasis):
        if shell_type0 < 0 or shell_type1 < 0:
            raise ValueError('A shell_type parameter can not be negative.')
        if max_nbasis < get_shell_nbasis(shell_type0):
            raise ValueError('max_nbasis to small for shell_type0.')
        if max_nbasis < get_shell_nbasis(shell_type1):
            raise ValueError('max_nbasis to small for shell_type1.')
        self._c_i2p.reset(shell_type0, shell_type1, max_nbasis)

    def inc(self):
        return self._c_i2p.inc()

    property fields:
        def __get__(self):
            return (
                self._c_i2p.n0[0], self._c_i2p.n0[1], self._c_i2p.n0[2],
                self._c_i2p.n1[0], self._c_i2p.n1[1], self._c_i2p.n1[2],
                self._c_i2p.offset
            )
