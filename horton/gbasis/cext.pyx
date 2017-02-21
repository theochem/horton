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

cimport libc.string
from libcpp.vector cimport vector
from cython.operator cimport dereference as deref

cimport boys
cimport cartpure
cimport common
cimport gbasis
cimport ints
cimport fns
cimport iter_gb
cimport iter_pow
cimport cholesky
cimport gbw

import atexit

from horton.log import log, biblio
from horton.matrix import LinalgFactory, CholeskyLinalgFactory
from horton.cext import compute_grid_nucpot
from horton.utils import typecheck_geo


__all__ = [
    # boys
    'boys_function', 'boys_function_array',
    # cartpure
    'cart_to_pure_low',
    #cholesky
    'compute_cholesky',
    # common
    'fac', 'fac2', 'binom', 'get_shell_nbasis', 'get_max_shell_type',
    'gpt_coeff', 'gb_overlap_int1d', 'nuclear_attraction_helper',
    'cit', 'jfac', 'dtaylor',
    # gbasis
    'gob_cart_normalization', 'gob_pure_normalization',
    'GOBasis',
    # gbw (testing)
    'get_2index_slice', 'compute_diagonal', 'select_2index',
    # ints
    'GB2OverlapIntegral', 'GB2KineticIntegral',
    'GB2ErfAttractionIntegral',
    'GB2GaussAttractionIntegral',
    'GB2NuclearAttractionIntegral','GB4Integral',
    'GB4ElectronRepulsionIntegralLibInt',
    'GB4ErfIntegralLibInt', 'GB4GaussIntegralLibInt',
    'GB4RAlphaIntegralLibInt',
    # fns
    'GB1DMGridDensityFn', 'GB1DMGridGradientFn', 'GB1DMGridGGAFn',
    'GB1DMGridKineticFn', 'GB1DMGridHessianFn', 'GB1DMGridMGGAFn',
    # iter_gb
    'IterGB1', 'IterGB2', 'IterGB4',
    # iter_pow
    'iter_pow1_inc', 'IterPow1', 'IterPow2',
]


#
# boys wrappers (for testing only)
#


def boys_function(long m, double t):
    return boys.boys_function(m, t)


def boys_function_array(long mmax, double t):
    cdef np.ndarray[double] output = np.zeros(mmax+1)
    boys.boys_function_array(mmax, t, &output[0])
    return output


#
# cartpure wrappers (for testing only)
#


def cart_to_pure_low(np.ndarray[double] work_cart not None,
                     np.ndarray[double] work_pure not None,
                     long shell_type, long nant, long npost):
    assert work_cart.flags['C_CONTIGUOUS']
    assert work_pure.flags['C_CONTIGUOUS']
    cartpure.cart_to_pure_low(
        &work_cart[0], &work_pure[0], shell_type, nant,
        npost
    )

#
# cholesky wrappers
#

def compute_cholesky(GOBasis gobasis, GB4Integral gb4int, double threshold=1e-8, lf = None):
    cdef gbw.GB4IntegralWrapper* gb4w = NULL
    cdef vector[np.float64_t]* vectors = NULL
    cdef np.npy_intp dims[3]
    cdef np.ndarray result

    cdef extern from "numpy/arrayobject.h":
        void PyArray_ENABLEFLAGS(np.ndarray arr, int flags)

    try:
        gb4w = new gbw.GB4IntegralWrapper(<gbasis.GOBasis*> gobasis._this,
                                          <ints.GB4Integral*> gb4int._this)
        vectors = new vector[np.float64_t]()
        nvec = cholesky.cholesky(gb4w, vectors, threshold)
        dims[0] = <np.npy_intp> nvec
        dims[1] = <np.npy_intp> gobasis.nbasis
        dims[2] = <np.npy_intp> gobasis.nbasis
        result = np.PyArray_SimpleNewFromData(3, dims, np.NPY_DOUBLE, &(deref(vectors)[0]))
        PyArray_ENABLEFLAGS(result, np.NPY_OWNDATA)
    finally:
        if gb4w is not NULL:
            del gb4w

    if lf is not None and isinstance(lf, CholeskyLinalgFactory):
        result_py = lf.create_four_index(gobasis.nbasis, array=result)
        return result_py
    else:
        return result

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


def gpt_coeff(long k, long n0, long n1, double pa, double pb):
    return common.gpt_coeff(k, n0, n1, pa, pb)


def gb_overlap_int1d(long n0, long n1, double pa, double pb, double inv_gamma):
    return common.gb_overlap_int1d(n0, n1, pa, pb, inv_gamma)


def nuclear_attraction_helper(np.ndarray[double, ndim=1] work_g not None,
                              long n0, long n1, double pa, double pb, double cp,
                              double gamma_inv):
    assert work_g.flags['C_CONTIGUOUS']
    assert work_g.shape[0] == n0+n1+1
    common.nuclear_attraction_helper(&work_g[0], n0, n1, pa, pb, cp, gamma_inv)

def cit(int i, double t, int m):
        return common.cit(i, t, m)

def jfac(int j, int n):
    return common.jfac(j, n)

def dtaylor(int n, double alpha, double t, double tfactor):
    return common.dtaylor(n, alpha, t, tfactor)


#
# gbasis wrappers
#


def gob_cart_normalization(double alpha, np.ndarray[long, ndim=1] n not None):
    assert n.flags['C_CONTIGUOUS']
    assert n.shape[0] == 3
    return gbasis.gob_cart_normalization(alpha, &n[0])


def gob_pure_normalization(double alpha, long l):
    return gbasis.gob_pure_normalization(alpha, l)


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

        self._centers.flags.writeable = True
        # Set arrays unwritable because:
        #   (i) derived properties will be stored
        #   (ii) consistency tests below are only performed once (below)
        # In principle, one can make the arrays writable again, but then one
        # is deliberately taking risks. A strict read-only buffer would be
        # ideal but is not possible. One can always pass a pointer to a
        # C library and start messing things up.
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

    def __init__(self, centers, shell_map, nprims, shell_types, alphas, con_coeffs):
        if self.__class__ == GBasis:
            raise NotImplementedError('GBasis is an abstract base class')
        self._log_init()

    def __dealloc__(self):
        del self._this

    @classmethod
    def concatenate(cls, *gbs):
        '''Concatenate multiple basis objects into a new one.

           **Arguments:** each argument is an instance of the same subclass of
           GBasis.
        '''
        # check if the classes match
        for gb in gbs:
            assert isinstance(gb, cls)

        # do the concatenation of each array properly
        centers = np.concatenate([gb.centers for gb in gbs])
        shell_map = []
        offset = 0
        for gb in gbs:
            shell_map.append(gb.shell_map + offset)
            offset += gb.ncenter
        shell_map = np.concatenate(shell_map)
        nprims = np.concatenate([gb.nprims for gb in gbs])
        shell_types = np.concatenate([gb.shell_types for gb in gbs])
        alphas = np.concatenate([gb.alphas for gb in gbs])
        con_coeffs = np.concatenate([gb.con_coeffs for gb in gbs])
        return cls(centers, shell_map, nprims, shell_types, alphas, con_coeffs)

    @classmethod
    def from_hdf5(cls, grp):
        return cls(
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

    # Array properties

    property centers:
        def __get__(self):
            return self._centers

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

    property ncenter:
        def __get__(self):
            return self.centers.shape[0]

    property nshell:
        def __get__(self):
            return self.shell_map.shape[0]

    property nprim_total:
        def __get__(self):
            return self.nprims.sum()

    # Other properties

    property nbasis:
        def __get__(self):
            return self._this.get_nbasis()

    property shell_lookup:
        def __get__(self):
            cdef np.npy_intp* shape = [self.nbasis]
            tmp = np.PyArray_SimpleNewFromData(1, shape,
                    np.NPY_LONG, <void*> self._this.get_shell_lookup())
            return tmp.copy()

    property basis_offsets:
        def __get__(self):
            cdef np.npy_intp* shape = [self.nshell]
            tmp = np.PyArray_SimpleNewFromData(1, shape, np.NPY_LONG,
                        <void*> self._this.get_basis_offsets())
            return tmp.copy()

    property nscales:
        def __get__(self):
            return self._this.get_nscales()

    property max_shell_type:
        def __get__(self):
            return self._this.get_max_shell_type()

    def _log_init(self):
        '''Write a summary of the basis to the screen logger'''
        if log.do_medium:
            log('Initialized: %s' % self)
            log.deflist([
                ('Number of basis functions', self.nbasis),
                ('Number of normalization constants', self.nscales),
                ('Maximum shell type', self.max_shell_type),
            ])
            shell_type_names = {
                0: 'S', 1: 'P', 2: 'Dc', 3: 'Fc', 4:'Gc', 5: 'Hc', 6: 'Ic',
                -2: 'Dp', -3: 'Fp', -4:'Gp', -5: 'Hp', -6: 'Ip',
            }
            descs = ['']*self.ncenter
            for i in xrange(self.nshell):
                icenter = self.shell_map[i]
                s = descs[icenter]
                name = shell_type_names[self.shell_types[i]]
                s += ' %s%i' % (name, self.nprims[i])
                descs[icenter] = s
            deflist = []
            for i in xrange(self.ncenter):
                deflist.append(('Center % 5i' % i, descs[i]))
            log.deflist(deflist)
            log.blank()

    def get_scales(self):
        # A **copy** of the scales is returned.
        cdef np.npy_intp shape[1]
        shape[0] = self.nscales
        tmp = np.PyArray_SimpleNewFromData(1, shape, np.NPY_DOUBLE, <void*> self._this.get_scales(0))
        return tmp.copy()

    # low-level compute routines, for debugging only
    def compute_grid_point1(self, np.ndarray[double, ndim=1] output not None,
                            np.ndarray[double, ndim=1] point not None,
                            GB1DMGridFn grid_fn not None):
        assert output.flags['C_CONTIGUOUS']
        assert output.shape[0] == self.nbasis
        assert point.flags['C_CONTIGUOUS']
        assert point.shape[0] == 3
        self._this.compute_grid_point1(&output[0], &point[0], grid_fn._this)

    def get_subset(self, ishells):
        '''Construct a sub basis set for a selection of shells

           **Argument:**

           ishells
                A list of indexes of shells to be retained in the sub basis set

           **Returns:** An instance of the same class as self containing only
           the basis functions of self that correspond to the select shells in
           the ``ishells`` list.
        '''
        # find the centers corresponding to ishells
        icenters = set([])
        for ishell in ishells:
            if ishell < 0:
                raise ValueError('ishell out of range: %i < 0' % ishell)
            if ishell >= self.nshell:
                raise ValueError('ishell out of range: %i >= %s' % (ishell, self.nshell))
            icenters.add(self.shell_map[ishell])
        icenters = sorted(icenters) # fix the order
        new_centers = self.centers[icenters]
        # construct the new shell_map, nprims, shell_types
        new_shell_map = np.zeros(len(ishells), int)
        new_nprims = np.zeros(len(ishells), int)
        new_shell_types = np.zeros(len(ishells), int)
        for new_ishell, ishell in enumerate(ishells):
            new_shell_map[new_ishell] = icenters.index(self.shell_map[ishell])
            new_nprims[new_ishell] = self.nprims[ishell]
            new_shell_types[new_ishell] = self.shell_types[ishell]
        # construct the new alphas and con_coeffs
        new_nprim_total = new_nprims.sum()
        new_alphas = np.zeros(new_nprim_total, float)
        new_con_coeffs = np.zeros(new_nprim_total, float)
        new_iprim = 0
        for new_ishell, ishell in enumerate(ishells):
            nprim = new_nprims[new_ishell]
            old_iprim = self.nprims[:ishell].sum()
            new_alphas[new_iprim:new_iprim+nprim] = self.alphas[old_iprim:old_iprim+nprim]
            new_con_coeffs[new_iprim:new_iprim+nprim] = self.con_coeffs[old_iprim:old_iprim+nprim]
            new_iprim += nprim
        # make a mapping between the indices of old and new basis functions
        ibasis_list = []
        for new_ishell, ishell in enumerate(ishells):
            ibasis_old = sum(common.get_shell_nbasis(self.shell_types[i]) for i in xrange(ishell))
            nbasis = common.get_shell_nbasis(self.shell_types[ishell])
            ibasis_list.extend(range(ibasis_old, ibasis_old+nbasis))
        ibasis_list = np.array(ibasis_list)
        # create the basis set object
        basis = self.__class__(new_centers, new_shell_map, new_nprims,
                               new_shell_types, new_alphas, new_con_coeffs)
        # return stuff
        return basis, ibasis_list

    def get_basis_atoms(self, coordinates):
        '''Return a list of atomic basis sets for a given geometry

           **Arguments:**

           coordinates
                An (N, 3) array with atomic coordinates, used to find the
                centers associated with atoms. An exact match of the Cartesian
                coordinates is required to properly select a shell.

           **Returns:** A list with one tuple for every atom: (gbasis,
           ibasis_list), where gbasis is a basis set object for the atom and
           ibasis_list is a list of basis set indexes that can be used to
           substitute results from the atomic basis set back into the molecular
           basis set. For example, when a density matrix for the atom is
           obtained and it needs to be plugged back into the molecular density
           matrix, one can do the following::

               mol_dm._array[ibasis_list,ibasis_list.reshape(-1,1)] = atom_dm._array
        '''
        result = []
        for c in coordinates:
            # find the corresponding center(s).
            icenters = []
            for icenter in xrange(self.ncenter):
                # require an exact match of the coordinates
                if (self.centers[icenter] == c).all():
                    icenters.append(icenter)
            icenters = set(icenters)
            # find the shells on these centers
            ishells = []
            for ishell in xrange(self.nshell):
                if self.shell_map[ishell] in icenters:
                    ishells.append(ishell)
            # construct a sub basis
            sub_basis, ibasis_list = self.get_subset(ishells)
            result.append((sub_basis, ibasis_list))
        return result


cdef class GOBasis(GBasis):
    def __cinit__(self, centers, shell_map, nprims, shell_types, alphas, con_coeffs):
        self._this = <gbasis.GBasis*> new gbasis.GOBasis(
            <double*>self._centers.data, <long*>self._shell_map.data,
            <long*>self._nprims.data, <long*>self._shell_types.data,
            <double*>self._alphas.data, <double*>self._con_coeffs.data,
            self._centers.shape[0], self._shell_types.shape[0], self._alphas.shape[0]
        )

    def check_matrix_coeffs(self, matrix, nocc=None):
        assert matrix.ndim == 2
        assert matrix.flags['C_CONTIGUOUS']
        assert matrix.shape[0] == self.nbasis
        assert matrix.shape[1] <= self.nbasis
        if nocc is not None:
            assert matrix.shape[1] >= nocc

    def check_matrix_two_index(self, matrix):
        assert matrix.ndim == 2
        assert matrix.flags['C_CONTIGUOUS']
        assert matrix.shape[0] == self.nbasis
        assert matrix.shape[1] == self.nbasis

    def check_matrix_four_index(self, matrix):
        assert matrix.ndim == 4
        assert matrix.flags['C_CONTIGUOUS']
        assert matrix.shape[0] == self.nbasis
        assert matrix.shape[1] == self.nbasis
        assert matrix.shape[2] == self.nbasis
        assert matrix.shape[3] == self.nbasis

    def compute_overlap(self, output):
        """Compute the overlap integrals in a Gaussian orbital basis

           **Arguments:**

           output
                When a ``TwoIndex`` instance is given, it is used as output
                argument and its contents are overwritten. When ``LinalgFactory``
                is given, it is used to construct the output ``TwoIndex``
                object. In both cases, the output two-index object is returned.

           **Returns:** ``TwoIndex`` object
        """
        # prepare the output array
        cdef np.ndarray[double, ndim=2] output_array
        if isinstance(output, LinalgFactory):
            lf = output
            output = lf.create_two_index(self.nbasis)
        output_array = output._array
        self.check_matrix_two_index(output_array)
        # call the low-level routine
        (<gbasis.GOBasis*>self._this).compute_overlap(&output_array[0, 0])
        # done
        return output

    def compute_kinetic(self, output):
        """Compute the kinetic energy integrals in a Gaussian orbital basis

           **Arguments:**

           output
                When a ``TwoIndex`` instance is given, it is used as output
                argument and its contents are overwritten. When ``LinalgFactory``
                is given, it is used to construct the output ``TwoIndex``
                object. In both cases, the output two-index object is returned.

           **Returns:** ``TwoIndex`` object
        """
        # prepare the output array
        cdef np.ndarray[double, ndim=2] output_array
        if isinstance(output, LinalgFactory):
            lf = output
            output = lf.create_two_index(self.nbasis)
        output_array = output._array
        self.check_matrix_two_index(output_array)
        # call the low-level routine
        (<gbasis.GOBasis*>self._this).compute_kinetic(&output_array[0, 0])
        # done
        return output

    def compute_nuclear_attraction(self,
                                   np.ndarray[double, ndim=2] coordinates not None,
                                   np.ndarray[double, ndim=1] charges not None,
                                   output):
        """Compute the nuclear attraction integral in a Gaussian orbital basis

           **Arguments:**

           coordinates
                A float array with shape (ncharge,3) with Cartesian coordinates
                of point charges that define the external field.

           charges
                A float array with shape (ncharge,) with the values of the
                charges.

           output
                When a ``TwoIndex`` instance is given, it is used as output
                argument and its contents are overwritten. When ``LinalgFactory``
                is given, it is used to construct the output ``TwoIndex``
                object. In both cases, the output two-index object is returned.

           **Returns:** ``TwoIndex`` object
        """
        # type checking
        assert coordinates.flags['C_CONTIGUOUS']
        assert charges.flags['C_CONTIGUOUS']
        ncharge, coordinates, charges = typecheck_geo(coordinates, None, charges, need_numbers=False)
        # prepare the output array
        cdef np.ndarray[double, ndim=2] output_array
        if isinstance(output, LinalgFactory):
            lf = output
            output = lf.create_two_index(self.nbasis)
        output_array = output._array
        self.check_matrix_two_index(output_array)
        # call the low-level routine
        (<gbasis.GOBasis*>self._this).compute_nuclear_attraction(
            &charges[0], &coordinates[0, 0], ncharge,
            &output_array[0, 0],
        )
        # done
        return output

    def compute_erf_attraction(self,
                                   np.ndarray[double, ndim=2] coordinates not None,
                                   np.ndarray[double, ndim=1] charges not None,
                                   output, double mu=0.0):
        r"""Compute the model nuclear attraction integral with the long-range potential

        The potential has the following form:

        .. math::
            v = \frac{\mathrm{erf}(\mu r)}{r}

        Parameters
        ----------
        coordinates: np.ndarray, shape = (ncharge,3)
            A float array with shape (ncharge,3) with Cartesian coordinates
            of point charges that define the external field.

        charges: np.ndarray, shape=(ncharge,)
             A float array with shape (ncharge,) with the values of the
             charges.

        output: TwoIndex
             When a ``TwoIndex`` instance is given, it is used as output
             argument and its contents are overwritten. When ``LinalgFactory``
             is given, it is used to construct the output ``TwoIndex``
             object. In both cases, the output two-index object is returned.

         mu : float
             Parameter for the erf(mu r)/r potential. Default is zero.

        Returns
        -------
        ``TwoIndex`` object

        Keywords: :index:`two-center integrals`
        """
        # type checking
        assert coordinates.flags['C_CONTIGUOUS']
        assert charges.flags['C_CONTIGUOUS']
        ncharge, coordinates, charges = typecheck_geo(coordinates, None, charges, need_numbers=False)
        # prepare the output array
        cdef np.ndarray[double, ndim=2] output_array
        if isinstance(output, LinalgFactory):
            lf = output
            output = lf.create_two_index(self.nbasis)
        output_array = output._array
        self.check_matrix_two_index(output_array)
        # call the low-level routine
        (<gbasis.GOBasis*>self._this).compute_erf_attraction(
            &charges[0], &coordinates[0, 0], ncharge,
            &output_array[0, 0], mu
        )
        # done
        return output

    def compute_gauss_attraction(self,
                                   np.ndarray[double, ndim=2] coordinates not None,
                                   np.ndarray[double, ndim=1] charges not None,
                                   output, double c=1.0, double alpha=1.0):
        r"""Compute the model nuclear attraction with a Gaussian potential.

        The potential has the following form:

        .. math::
            v = c \exp(-\alpha r^2)

        Parameters
        ----------
        output : TwoIndex
            When a ``DenseTwoIndex`` object is given, it is used as output argument and
            its contents are overwritten.
        c : float
            Coefficient of the gaussian.
        alpha : float
            Exponential parameter of the gaussian.

        Returns
        -------
        output

        Keywords: :index:`two-center integrals`
        """
        # type checking
        assert coordinates.flags['C_CONTIGUOUS']
        assert charges.flags['C_CONTIGUOUS']
        ncharge, coordinates, charges = typecheck_geo(coordinates, None, charges, need_numbers=False)
        # prepare the output array
        cdef np.ndarray[double, ndim=2] output_array
        if isinstance(output, LinalgFactory):
            lf = output
            output = lf.create_two_index(self.nbasis)
        output_array = output._array
        self.check_matrix_two_index(output_array)
        # call the low-level routine
        (<gbasis.GOBasis*>self._this).compute_gauss_attraction(
            &charges[0], &coordinates[0, 0], ncharge,
            &output_array[0, 0], c, alpha
        )
        # done
        return output

    def compute_multipole_moment(self,
                    np.ndarray[long, ndim=1] xyz not None,
                    np.ndarray[double, ndim=1] center not None,
                    output):
        """Compute the (multipole) moment integrals in a Gaussian orbital basis.

        Calculates the integral < gto_a | (x - C_x)^l (y - C_y)^m (z - C_z)^n | gto_b >

        Parameters
        ----------
        xyz : numpy-array of int, shape=(3,).
            A integer (long) numpy-array with shape (3,) with the powers of x,y,z in the
            integrals.
        center : np.ndarray, shape = (3,)
            A numpy array of shape (3,) with the center [C_x, C_y, C_z] around which the
            moment integral is computed.
        output : ``TwoIndex`` or ``LinalgFactory`` object
            When a ``TwoIndex`` instance is given, it is used as output argument and its
            contents are overwritten. When ``LinalgFactory`` is given, it is used to
            construct the output ``TwoIndex`` object. In both cases, the output two-index
            object is returned.

        Returns
        -------
        output : ``TwoIndex`` object
            The values of the integrals.
        """
        # type checking
        assert xyz.flags['C_CONTIGUOUS']
        assert xyz.min() >= 0
        assert xyz.sum() >= 0
        assert center.flags['C_CONTIGUOUS']
        # prepare the output array
        cdef np.ndarray[double, ndim=2] output_array
        if isinstance(output, LinalgFactory):
            lf = output
            output = lf.create_two_index(self.nbasis)
        output_array = output._array
        self.check_matrix_two_index(output_array)
        # call the low-level routine
        (<gbasis.GOBasis*>self._this).compute_multipole_moment(
            &xyz[0], &center[0], &output_array[0, 0])
        # done
        return output

    def compute_electron_repulsion(self, output):
        r'''Compute electron-electron repulsion integrals

        The potential has the following form:

        .. math::
            v = \frac{1}{r}

        Parameters
        ----------
        output : FourIndex
            When a ``DenseFourIndex`` object is given, it is used as output argument and
            its contents are overwritten. When a ``DenseLinalgFactory`` or
            ``CholeskyLinalgFactory`` is given, it is used to construct the four-index
            object in which the integrals are stored.

        Returns
        -------
        output

        Keywords: :index:`ERI`, :index:`four-center integrals`
        '''
        biblio.cite('valeev2014',
                 'the efficient implementation of four-center electron repulsion integrals')
        if isinstance(output, CholeskyLinalgFactory):
            lf = output
            output = compute_cholesky(self, GB4ElectronRepulsionIntegralLibInt(self.max_shell_type), lf=lf)
            return output
        # prepare the output array
        cdef np.ndarray[double, ndim=4] output_array
        if isinstance(output, LinalgFactory):
            lf = output
            output = lf.create_four_index(self.nbasis)
        output_array = output._array
        self.check_matrix_four_index(output_array)
        # call the low-level routine
        (<gbasis.GOBasis*>self._this).compute_electron_repulsion(&output_array[0, 0, 0, 0])
        # done
        return output

    def compute_erf_repulsion(self, output, double mu=0.0):
        r"""Compute short-range electron repulsion integrals.

        The potential has the following form:

        .. math::
            v = \frac{\mathrm{erf}(\mu r)}{r}

        Parameters
        ----------
        output : FourIndex
            When a ``DenseFourIndex`` object is given, it is used as output argument and
            its contents are overwritten. When a ``DenseLinalgFactory`` or
            ``CholeskyLinalgFactory`` is given, it is used to construct the four-index
            object in which the integrals are stored.
        mu : float
            Parameter for the erf(mu r)/r potential. Default is zero.

        Returns
        -------
        output

        Keywords: :index:`ERI`, :index:`four-center integrals`
        """
        biblio.cite('valeev2014',
                 'the efficient implementation of four-center electron repulsion integrals')
        biblio.cite('ahlrichs2006',
                 'the methodology to implement various types of four-center integrals.')
        if isinstance(output, CholeskyLinalgFactory):
            lf = output
            output = compute_cholesky(self, GB4ErfIntegralLibInt(self.max_shell_type, mu), lf=lf)
            return output
        # prepare the output array
        cdef np.ndarray[double, ndim=4] output_array
        if isinstance(output, LinalgFactory):
            lf = output
            output = lf.create_four_index(self.nbasis)
        output_array = output._array
        self.check_matrix_four_index(output_array)
        # call the low-level routine
        (<gbasis.GOBasis*>self._this).compute_erf_repulsion(&output_array[0, 0, 0, 0], mu)
        # done
        return output

    def compute_gauss_repulsion(self, output, double c=1.0, double alpha=1.0):
        r"""Compute gaussian repulsion four-center integrals.

        The potential has the following form:

        .. math::
            v = c \exp(-\alpha r^2)

        Parameters
        ----------
        output : FourIndex
            When a ``DenseFourIndex`` object is given, it is used as output argument and
            its contents are overwritten. When a ``DenseLinalgFactory`` or
            ``CholeskyLinalgFactory`` is given, it is used to construct the four-index
            object in which the integrals are stored.
        c : float
            Coefficient of the gaussian.
        alpha : float
            Exponential parameter of the gaussian.

        Returns
        -------
        output

        Keywords: :index:`ERI`, :index:`four-center integrals`
        """
        biblio.cite('valeev2014',
                 'the efficient implementation of four-center electron repulsion integrals')
        biblio.cite('ahlrichs2006',
                 'the methodology to implement various types of four-center integrals.')
        biblio.cite('gill1996',
                 'four-center integrals with a Gaussian interaction potential.')
        biblio.cite('toulouse2004',
                 'four-center integrals with a Gaussian interaction potential.')
        if isinstance(output, CholeskyLinalgFactory):
            lf = output
            output = compute_cholesky(self, GB4GaussIntegralLibInt(self.max_shell_type, c, alpha), lf=lf)
            return output
        # prepare the output array
        cdef np.ndarray[double, ndim=4] output_array
        if isinstance(output, LinalgFactory):
            lf = output
            output = lf.create_four_index(self.nbasis)
        output_array = output._array
        self.check_matrix_four_index(output_array)
        # call the low-level routine
        (<gbasis.GOBasis*>self._this).compute_gauss_repulsion(&output_array[0, 0, 0, 0], c, alpha)
        # done
        return output

    def compute_ralpha_repulsion(self, output, double alpha=-1.0):
        r"""Compute r^alpha repulsion four-center integrals.

        The potential has the following form:

        .. math::
            v = r^{\alpha}

        with :math:`\alpha > -3`.

        Parameters
        ----------

        output : FourIndex
            When a ``DenseFourIndex`` object is given, it is used as output argument and
            its contents are overwritten. When a ``DenseLinalgFactory`` or
            ``CholeskyLinalgFactory`` is given, it is used to construct the four-index
            object in which the integrals are stored.
        alpha : float
            The power of r in the interation potential.

        Returns
        -------
        output

        Keywords: :index:`ERI`, :index:`four-center integrals`
        """
        biblio.cite('valeev2014',
                 'the efficient implementation of four-center electron repulsion integrals')
        biblio.cite('ahlrichs2006',
                 'the methodology to implement various types of four-center integrals.')
        if isinstance(output, CholeskyLinalgFactory):
            lf = output
            output = compute_cholesky(self, GB4RAlphaIntegralLibInt(self.max_shell_type, alpha), lf=lf)
            return output
        cdef np.ndarray[double, ndim=4] output_array
        if isinstance(output, LinalgFactory):
            lf = output
            output = lf.create_four_index(self.nbasis)
        output_array = output._array
        self.check_matrix_four_index(output_array)
        # call the low-level routine
        (<gbasis.GOBasis*>self._this).compute_ralpha_repulsion(&output_array[0, 0, 0, 0], alpha)
        # done
        return output

    def compute_grid_orbitals_exp(self, exp,
                                  np.ndarray[double, ndim=2] points not None,
                                  np.ndarray[long, ndim=1] iorbs not None,
                                  np.ndarray[double, ndim=2] output=None):
        """Compute the orbitals on a grid for a given set of expansion coefficients.

        **Warning:** the results are added to the output array!

        Parameters
        ----------

        exp : DenseExpansion
            Orbitals.
        points : np.ndarray, shape=(npoint, 3), dtype=float
            Cartesian grid points.
        iorbs : np.ndarray, shape=(n,), dtype=int
            Indexes of the orbitals to be computed. When not given, the orbitals with a
            non-zero occupation number are computed.
        output : np.ndarray, shape=(npoint, n), dtype=float
            An output array. The results are added to this array. When not given, an
            output array is allocated.

        Returns
        -------
        output : np.ndarray, shape=(npoint, n), dtype=float
            the output array. (It is allocated when not given.)
        """
        # Do some type checking
        cdef np.ndarray[double, ndim=2] coeffs = exp.coeffs
        self.check_matrix_coeffs(coeffs)
        nfn = coeffs.shape[1]
        assert points.flags['C_CONTIGUOUS']
        npoint = points.shape[0]
        assert points.shape[1] == 3
        assert iorbs.flags['C_CONTIGUOUS']
        norb = iorbs.shape[0]
        if output is None:
            output = np.zeros((npoint, norb), float)
        else:
            assert output.flags['C_CONTIGUOUS']
            assert output.shape[0] == npoint
            assert output.shape[1] == norb
        # compute
        (<gbasis.GOBasis*>self._this).compute_grid1_exp(
            nfn, &coeffs[0, 0], npoint, &points[0, 0],
            norb, &iorbs[0], &output[0, 0])
        return output

    def _compute_grid1_dm(self, dm, np.ndarray[double, ndim=2] points not None,
                          GB1DMGridFn grid_fn not None, np.ndarray output not None,
                          double epsilon=0):
        """Compute some density function on a grid for a given density matrix.

        **Warning:** the results are added to the output array! This may be useful to
        combine results from different spin components.

        Parameters
        ----------

        dm : DenseTwoIndex
            Density matrix, assumed to be symmetric.
        points : np.ndarray, shape=(npoint, 3), dtype=float
            Cartesian grid points.
        grid_fn : GB1DMGridFn
            Implements the function to be evaluated on the grid.
        output : np.ndarray, shape=(npoint, n), dtype=float
            Output array. The second dimension depends on grid_fn.
        epsilon : float
            Allow errors on the density of this magnitude for the sake of
            efficiency. Some grid_fn implementations may ignore this.
        """
        # Get the array of the density matrix
        cdef np.ndarray[double, ndim=2] dmar = dm._array
        self.check_matrix_two_index(dmar)

        # Get the maximum of the absolute value over the rows
        cdef np.ndarray[double, ndim=1] dmmaxrow = np.abs(dmar).max(axis=0)

        # Check the output array
        assert output.flags['C_CONTIGUOUS']
        assert output.dtype == np.double
        npoint = output.shape[0]
        if grid_fn.dim_output == 1:
            assert output.ndim == 1
        else:
            assert output.ndim == 2
            assert output.shape[1] == grid_fn.dim_output

        # Check the points array
        assert points.flags['C_CONTIGUOUS']
        assert points.shape[0] == npoint
        assert points.shape[1] == 3

        # Go!
        (<gbasis.GOBasis*>self._this).compute_grid1_dm(
            &dmar[0, 0], npoint, &points[0, 0],
            grid_fn._this, <double*>np.PyArray_DATA(output), epsilon,
            &dmmaxrow[0])

    def compute_grid_density_dm(self, dm,
                                np.ndarray[double, ndim=2] points not None,
                                np.ndarray[double, ndim=1] output=None,
                                double epsilon=0):
        """Compute the electron density on a grid for a given density matrix.

        **Warning:** the results are added to the output array! This may be useful to
        combine results from different spin components.

        Parameters
        ----------

        dm : DenseTwoIndex
            Density matrix, assumed to be symmetric.
        points : np.ndarray, shape=(npoint, 3), dtype=float
            Cartesian grid points.
        output : np.ndarray, shape=(npoint,), dtype=float
            Output array. When not given, it is allocated and returned.
        epsilon : float
            Allow errors on the density of this magnitude for the sake of
            efficiency. Some grid_fn implementations may ignore this.

        Returns
        -------
        output : np.ndarray, shape=(npoint,), dtype=float
            The output array.
        """
        if output is None:
            output = np.zeros(points.shape[0])
        self._compute_grid1_dm(dm, points, GB1DMGridDensityFn(self.max_shell_type), output, epsilon)
        return output

    def compute_grid_gradient_dm(self, dm,
                                 np.ndarray[double, ndim=2] points not None,
                                 np.ndarray[double, ndim=2] output=None):
        """Compute the electron density gradient on a grid for a given density matrix.

        **Warning:** the results are added to the output array! This may be useful to
        combine results from different spin components.

        Parameters
        ----------

        dm : DenseTwoIndex
            Density matrix, assumed to be symmetric.
        points : np.ndarray, shape=(npoint, 3), dtype=float
            Cartesian grid points.
        output : np.ndarray, shape=(npoint, 3), dtype=float
            Output array. When not given, it is allocated and returned.

        Returns
        -------
        output : np.ndarray, shape=(npoint, 3), dtype=float
            The output array.
        """
        if output is None:
            output = np.zeros((points.shape[0], 3), float)
        self._compute_grid1_dm(dm, points, GB1DMGridGradientFn(self.max_shell_type), output)
        return output

    def compute_grid_gga_dm(self, dm,
                            np.ndarray[double, ndim=2] points not None,
                            np.ndarray[double, ndim=2] output=None):
        """Compute the electron density and gradient on a grid for a given density matrix.

        **Warning:** the results are added to the output array! This may be useful to
        combine results from different spin components.

        Parameters
        ----------

        dm : DenseTwoIndex
            Density matrix, assumed to be symmetric.
        points : np.ndarray, shape=(npoint, 3), dtype=float
            Cartesian grid points.
        output : np.ndarray, shape=(npoint, 4), dtype=float
            Output array. When not given, it is allocated and returned. The first column
            contains the density. The last three columns contain the gradient.

        Returns
        -------
        output : np.ndarray, shape=(npoint, 4), dtype=float
            The output array.
        """
        if output is None:
            output = np.zeros((points.shape[0], 4), float)
        self._compute_grid1_dm(dm, points, GB1DMGridGGAFn(self.max_shell_type), output)
        return output

    def compute_grid_kinetic_dm(self, dm,
                                np.ndarray[double, ndim=2] points not None,
                                np.ndarray[double, ndim=1] output=None):
        """Compute the kinetic energy density on a grid for a given density matrix.

        **Warning:** the results are added to the output array! This may be useful to
        combine results from different spin components.

        Parameters
        ----------

        dm : DenseTwoIndex
            Density matrix, assumed to be symmetric.
        points : np.ndarray, shape=(npoint, 3), dtype=float
            Cartesian grid points.
        output : np.ndarray, shape=(npoint,), dtype=float
            Output array. When not given, it is allocated and returned.

        Returns
        -------
        output : np.ndarray, shape=(npoint,), dtype=float
            The output array.
        """
        if output is None:
            output = np.zeros((points.shape[0],), float)
        self._compute_grid1_dm(dm, points, GB1DMGridKineticFn(self.max_shell_type), output)
        return output

    def compute_grid_hessian_dm(self, dm,
                                np.ndarray[double, ndim=2] points not None,
                                np.ndarray[double, ndim=2] output=None):
        """Compute the electron density Hessian on a grid for a given density matrix.

        **Warning:** the results are added to the output array! This may be useful to
        combine results from different spin components.

        Parameters
        ----------

        dm : DenseTwoIndex
            Density matrix, assumed to be symmetric.
        points : np.ndarray, shape=(npoint, 3), dtype=float
            Cartesian grid points.
        output : np.ndarray, shape=(npoint, 6), dtype=float
            Output array. When not given, it is allocated and returned. The columns are
            assigned as follows:

            * 0: element (0, 0) of the Hessian
            * 1: element (0, 1) of the Hessian
            * 2: element (0, 2) of the Hessian
            * 3: element (1, 1) of the Hessian
            * 4: element (1, 2) of the Hessian
            * 5: element (2, 2) of the Hessian

        Returns
        -------
        output : np.ndarray, shape=(npoint, 6), dtype=float
            The output array.
        """
        if output is None:
            output = np.zeros((points.shape[0], 6), float)
        self._compute_grid1_dm(dm, points, GB1DMGridHessianFn(self.max_shell_type), output)
        return output

    def compute_grid_mgga_dm(self, dm,
                             np.ndarray[double, ndim=2] points not None,
                             np.ndarray[double, ndim=2] output=None):
        """Compute the MGGA quantities for a given density matrix.

        **Warning:** the results are added to the output array! This may be useful to
        combine results from different spin components.

        This includes the density, the gradient, the Laplacian and the
        kinetic energy density.

        Parameters
        ----------

        dm : DenseTwoIndex
            Density matrix, assumed to be symmetric.
        points : np.ndarray, shape=(npoint, 3), dtype=float
            Cartesian grid points.
        output : np.ndarray, shape=(npoint, 6), dtype=float
            Output array. When not given, it is allocated and returned. The assignment of
            the columns is as follows:

            * 0: density
            * 1: gradient x
            * 2: gradient y
            * 3: gradient z
            * 4: laplacian
            * 5: kinetic energy density


        Returns
        -------
        output : np.ndarray, shape=(npoint, 6), dtype=float
            The output array.
        """
        if output is None:
            output = np.zeros((points.shape[0], 6), float)
        self._compute_grid1_dm(dm, points, GB1DMGridMGGAFn(self.max_shell_type), output)
        return output

    def compute_grid_hartree_dm(self, dm,
                                np.ndarray[double, ndim=2] points not None,
                                np.ndarray[double, ndim=1] output=None):
        """Compute the Hartree potential on a grid for a given density matrix.

        **Warning:** the results are added to the output array! This may be useful to
        combine results from different spin components.

        Parameters
        ----------

        dm : DenseTwoIndex
            Density matrix, assumed to be symmetric.
        points : np.ndarray, shape=(npoint, 3), dtype=float
            Cartesian grid points.
        output : np.ndarray, shape=(npoint,), dtype=float
            Output array. When not given, it is allocated and returned.


        Returns
        -------
        output : np.ndarray, shape=(npoint,), dtype=float
            The output array.
        """
        # type checking
        cdef np.ndarray[double, ndim=2] dmar = dm._array
        self.check_matrix_two_index(dmar)
        assert points.flags['C_CONTIGUOUS']
        npoint = points.shape[0]
        assert points.shape[1] == 3
        if output is None:
            output = np.zeros(npoint, float)
        else:
            assert output.flags['C_CONTIGUOUS']
            assert output.shape[0] == npoint
        # compute
        (<gbasis.GOBasis*>self._this).compute_grid2_dm(
            &dmar[0, 0], npoint, &points[0, 0],
            &output[0])
        return output

    def compute_grid_esp_dm(self, dm,
                            np.ndarray[double, ndim=2] coordinates not None,
                            np.ndarray[double, ndim=1] charges not None,
                            np.ndarray[double, ndim=2] points not None,
                            np.ndarray[double, ndim=1] output=None):
        """Compute the electrostatic potential on a grid for a given density matrix.

        **Warning:** the results are added to the output array! This may be useful to
        combine results from different spin components.

        Parameters
        ----------

        dm : DenseTwoIndex
            Density matrix, assumed to be symmetric.
        coordinates : np.ndarray, shape=(natom, 3), dtype=float
            Cartesian coordinates of the atoms.
        charges : np.ndarray, shape=(natom,), dtype=float
            Atomic charges.
        points : np.ndarray, shape=(npoint, 3), dtype=float
            Cartesian grid points.
        output : np.ndarray, shape=(npoint,), dtype=float
            Output array. When not given, it is allocated and returned.


        Returns
        -------
        output : np.ndarray, shape=(npoint,), dtype=float
            The output array.
        """
        output = self.compute_grid_hartree_dm(dm, points, output)
        output *= -1
        compute_grid_nucpot(coordinates, charges, points, output)
        return output

    def _compute_grid1_fock(self, np.ndarray[double, ndim=2] points not None,
                           np.ndarray[double, ndim=1] weights not None,
                           np.ndarray pots not None,
                           GB1DMGridFn grid_fn not None, fock):
        """Compute a Fock operator from a some sort of potential.

        **Warning:** the results are added to the Fock operator!

        Parameters
        ----------

        points : np.ndarray, shape=(npoint, 3), dtype=float
            Cartesian grid points.
        weights : np.ndarray, shape=(npoint,), dtype=float
            Integration weights.
        pots : np.ndarray, shape=(npoint, n)
            Derivative of the energy toward the density-related quantities
            at all grid points. The number of columns depends on grid_fn.
        grid_fn : GB1DMGridFn
            Implements the function to be evaluated on the grid.
        fock : DenseTwoIndex
            Output dense two-index object.
        """
        cdef np.ndarray[double, ndim=2] output = fock._array
        self.check_matrix_two_index(output)
        assert points.flags['C_CONTIGUOUS']
        npoint = points.shape[0]
        assert points.shape[1] == 3
        assert weights.flags['C_CONTIGUOUS']
        assert npoint == weights.shape[0]
        assert pots.strides[0] % 8 == 0
        pot_stride = pots.strides[0]/8
        assert npoint == pots.shape[0]
        if grid_fn.dim_output == 1:
            assert pots.ndim == 1
        else:
            assert pots.ndim == 2
            assert pots.shape[1] == grid_fn.dim_output
            assert pots.strides[1] % 8 == 0
            pot_stride *= (pots.strides[1] / 8)
        (<gbasis.GOBasis*>self._this).compute_grid1_fock(
            npoint, &points[0, 0], &weights[0],
            pot_stride, <double*>np.PyArray_DATA(pots),
            grid_fn._this, &output[0, 0])

    def compute_grid_density_fock(self, np.ndarray[double, ndim=2] points not None,
                                  np.ndarray[double, ndim=1] weights not None,
                                  np.ndarray[double, ndim=1] pots not None, fock):
        """Compute a Fock operator from a density potential.

        **Warning:** the results are added to the Fock operator!

        Parameters
        ----------

        points : np.ndarray, shape=(npoint, 3), dtype=float
            Cartesian grid points.
        weights : np.ndarray, shape=(npoint,), dtype=float
            Integration weights.
        pots : np.ndarray, shape=(npoint,), dtype=float
            Derivative of the energy toward the density at all grid points.
        fock : DenseTwoIndex
            Output dense two-index object.
        """
        self._compute_grid1_fock(points, weights, pots, GB1DMGridDensityFn(self.max_shell_type), fock)

    def compute_grid_gradient_fock(self, np.ndarray[double, ndim=2] points not None,
                                   np.ndarray[double, ndim=1] weights not None,
                                   np.ndarray[double, ndim=2] pots not None, fock):
        """Compute a Fock operator from a density gradient potential.

        **Warning:** the results are added to the Fock operator!

        Parameters
        ----------

        points : np.ndarray, shape=(npoint, 3), dtype=float
            Cartesian grid points.
        weights : np.ndarray, shape=(npoint,), dtype=float
            Integration weights.
        pots : np.ndarray, shape=(npoint, 3), dtype=float
            Derivative of the energy toward the density gradient components at all grid
            points.
        fock : DenseTwoIndex
            Output dense two-index object.
        """
        self._compute_grid1_fock(points, weights, pots, GB1DMGridGradientFn(self.max_shell_type), fock)

    def compute_grid_gga_fock(self, np.ndarray[double, ndim=2] points not None,
                                    np.ndarray[double, ndim=1] weights not None,
                                    np.ndarray[double, ndim=2] pots not None, fock):
        """Compute a Fock operator from GGA potential data.

        **Warning:** the results are added to the Fock operator!

        Parameters
        ----------

        points : np.ndarray, shape=(npoint, 3), dtype=float
            Cartesian grid points.
        weights : np.ndarray, shape=(npoint,), dtype=float
            Integration weights.
        pots : np.ndarray, shape=(npoint, 4), dtype=float
            Derivative of the energy toward GGA ingredients (density and gradient) at all
            grid points.
        fock : DenseTwoIndex
            Output dense two-index object.
        """
        # To be replaced by something more efficient
        self._compute_grid1_fock(points, weights, pots, GB1DMGridGGAFn(self.max_shell_type), fock)

    def compute_grid_kinetic_fock(self, np.ndarray[double, ndim=2] points not None,
                                  np.ndarray[double, ndim=1] weights not None,
                                  np.ndarray[double, ndim=1] pots not None, fock):
        """Compute a Fock operator from a kientic-energy-density potential.

        **Warning:** the results are added to the Fock operator!

        Parameters
        ----------

        points : np.ndarray, shape=(npoint, 3), dtype=float
            Cartesian grid points.
        weights : np.ndarray, shape=(npoint,), dtype=float
            Integration weights.
        pots : np.ndarray, shape=(npoint,), dtype=float
            Derivative of the energy toward the kinetic energy density at all grid points.
        fock : DenseTwoIndex
            Output dense two-index object.
        """
        self._compute_grid1_fock(points, weights, pots, GB1DMGridKineticFn(self.max_shell_type), fock)

    def compute_grid_hessian_fock(self, np.ndarray[double, ndim=2] points not None,
                                  np.ndarray[double, ndim=1] weights not None,
                                  np.ndarray[double, ndim=2] pots not None, fock):
        """Compute a Fock operator from a density hessian potential.

        **Warning:** the results are added to the Fock operator!

        Parameters
        ----------

        points : np.ndarray, shape=(npoint, 3), dtype=float
            Cartesian grid points.
        weights : np.ndarray, shape=(npoint,), dtype=float
            Integration weights.
        pots : np.ndarray, shape=(npoint, 6), dtype=float
            Derivative of energy toward components of the Hessian at all grid points. The
            columns are assigned as follows:

            * 0: element (0, 0) of the Hessian
            * 1: element (0, 1) of the Hessian
            * 2: element (0, 2) of the Hessian
            * 3: element (1, 1) of the Hessian
            * 4: element (1, 2) of the Hessian
            * 5: element (2, 2) of the Hessian
        fock : DenseTwoIndex
            Output dense two-index object.
        """
        self._compute_grid1_fock(points, weights, pots, GB1DMGridHessianFn(self.max_shell_type), fock)

    def compute_grid_mgga_fock(self, np.ndarray[double, ndim=2] points not None,
                                np.ndarray[double, ndim=1] weights not None,
                                np.ndarray[double, ndim=2] pots not None, fock):
        """Compute a Fock operator from MGGA potential data.

        **Warning:** the results are added to the Fock operator!

        Parameters
        ----------

        points : np.ndarray, shape=(npoint, 3), dtype=float
            Cartesian grid points.
        weights : np.ndarray, shape=(npoint,), dtype=float
            Integration weights.
        pots : np.ndarray, shape=(npoint, 6), dtype=float
            Derivative of the energy toward density, gradient, Laplacian and kinetic
            energy density. The assignment of the columns is as follows:

            * 0: density
            * 1: gradient x
            * 2: gradient y
            * 3: gradient z
            * 4: laplacian
            * 5: kinetic energy density
        fock : DenseTwoIndex
            Output dense two-index object.
        """
        self._compute_grid1_fock(points, weights, pots, GB1DMGridMGGAFn(self.max_shell_type), fock)


#
# gbw wrappers
#

def select_2index(GOBasis gobasis, long index0, long index2):
    assert index0 >= 0 and index0 < gobasis.nbasis
    assert index2 >= 0 and index2 < gobasis.nbasis

    cdef ints.GB4ElectronRepulsionIntegralLibInt* gb4int = NULL
    cdef gbw.GB4IntegralWrapper* gb4w = NULL

    cdef long pbegin0
    cdef long pend0
    cdef long pbegin2
    cdef long pend2

    try:
        gb4int = new ints.GB4ElectronRepulsionIntegralLibInt(
                            gobasis.max_shell_type)
        gb4w = new gbw.GB4IntegralWrapper((<gbasis.GOBasis* > gobasis._this),
                            <ints.GB4Integral*> gb4int)
        gb4w.select_2index(index0, index2, &pbegin0, &pend0, &pbegin2, &pend2)
    finally:
        if gb4int is not NULL:
            del gb4int
        if gb4w is not NULL:
            del gb4w
    return pbegin0, pend0, pbegin2, pend2


def compute_diagonal(GOBasis gobasis, np.ndarray[double, ndim=2] diagonal not
        None):
    cdef ints.GB4ElectronRepulsionIntegralLibInt* gb4int = NULL
    cdef gbw.GB4IntegralWrapper* gb4w = NULL
    cdef np.ndarray[double, ndim=2] output
    output = diagonal

    try:
        gb4int = new ints.GB4ElectronRepulsionIntegralLibInt(
                            gobasis.max_shell_type)
        gb4w = new gbw.GB4IntegralWrapper((<gbasis.GOBasis* > gobasis._this),
                            <ints.GB4Integral*> gb4int)
        gb4w.compute_diagonal(&output[0, 0])

    finally:
        if gb4int is not NULL:
            del gb4int
        if gb4w is not NULL:
            del gb4w

def get_2index_slice(GOBasis gobasis, long index0, long index2,
                        np.ndarray[double, ndim=2] slice not None):
    cdef ints.GB4ElectronRepulsionIntegralLibInt* gb4int = NULL
    cdef gbw.GB4IntegralWrapper* gb4w = NULL
    assert slice.flags['C_CONTIGUOUS']
    assert slice.shape[0] == gobasis.nbasis
    assert slice.shape[1] == gobasis.nbasis

    cdef long pbegin0
    cdef long pend0
    cdef long pbegin2
    cdef long pend2
    cdef double* output
    try:
        gb4int = new ints.GB4ElectronRepulsionIntegralLibInt(
                            gobasis.max_shell_type)
        gb4w = new gbw.GB4IntegralWrapper((<gbasis.GOBasis* > gobasis._this),
                            <ints.GB4Integral*> gb4int)
        gb4w.select_2index(index0, index2, &pbegin0, &pend0, &pbegin2, &pend2)
        gb4w.compute()
        output = gb4w.get_2index_slice(index0, index2)
        print output[0]
        print sizeof(double)*gobasis.nbasis*gobasis.nbasis
        libc.string.memcpy(&slice[0,0], output,
                sizeof(double)*gobasis.nbasis*gobasis.nbasis)
        print slice[0,0]

    finally:
        if gb4int is not NULL:
            del gb4int
        if gb4w is not NULL:
            del gb4w


#
# ints wrappers (for testing only)
#


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
              np.ndarray[double, ndim=1] r0 not None,
              np.ndarray[double, ndim=1] r1 not None):
        assert r0.flags['C_CONTIGUOUS']
        assert r0.shape[0] == 3
        assert r1.flags['C_CONTIGUOUS']
        assert r1.shape[0] == 3
        self._this.reset(shell_type0, shell_type1, <double*>r0.data, <double*>r1.data)

    def add(self, double coeff, double alpha0, double alpha1,
            np.ndarray[double, ndim=1] scales0 not None,
            np.ndarray[double, ndim=1] scales1 not None):
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
    '''Wrapper for ints.GB2KineticIntegral, for testing only'''

    def __cinit__(self, long max_nbasis):
        self._this = <ints.GB2Integral*>(new ints.GB2KineticIntegral(max_nbasis))


cdef class GB2NuclearAttractionIntegral(GB2Integral):
    '''Wrapper for ints.GB2NuclearAttractionIntegral, for testing only'''
    # make an additional reference to these arguments to avoid deallocation
    cdef np.ndarray _charges
    cdef np.ndarray _centers

    def __cinit__(self, long max_nbasis,
                  np.ndarray[double, ndim=1] charges not None,
                  np.ndarray[double, ndim=2] centers not None):
        assert charges.flags['C_CONTIGUOUS']
        cdef long ncharge = charges.shape[0]
        assert centers.flags['C_CONTIGUOUS']
        assert centers.shape[0] == ncharge
        self._charges = charges
        self._centers = centers
        self._this = <ints.GB2Integral*>(new ints.GB2NuclearAttractionIntegral(
           max_nbasis, &charges[0], &centers[0, 0], ncharge
        ))


cdef class GB2ErfAttractionIntegral(GB2Integral):
    '''Wrapper for ints.GB2ErfAttractionIntegral, for testing only'''
    # make an additional reference to these arguments to avoid deallocation
    cdef np.ndarray _charges
    cdef np.ndarray _centers

    def __cinit__(self, long max_nbasis,
                  np.ndarray[double, ndim=1] charges not None,
                  np.ndarray[double, ndim=2] centers not None, double mu):
        assert charges.flags['C_CONTIGUOUS']
        cdef long ncharge = charges.shape[0]
        assert centers.flags['C_CONTIGUOUS']
        assert centers.shape[0] == ncharge
        self._charges = charges
        self._centers = centers
        self._this = <ints.GB2Integral*>(new ints.GB2ErfAttractionIntegral(
            max_nbasis, &charges[0], &centers[0, 0], ncharge, mu
        ))

    property mu:
        def __get__(self):
            return (<ints.GB2ErfAttractionIntegral*>self._this).get_mu()


cdef class GB2GaussAttractionIntegral(GB2Integral):
    '''Wrapper for ints.GB2GaussAttractionIntegral, for testing only'''
    # make an additional reference to these arguments to avoid deallocation
    cdef np.ndarray _charges
    cdef np.ndarray _centers

    def __cinit__(self, long max_nbasis,
                  np.ndarray[double, ndim=1] charges not None,
                  np.ndarray[double, ndim=2] centers not None, double c,
                  double alpha):
        assert charges.flags['C_CONTIGUOUS']
        cdef long ncharge = charges.shape[0]
        assert centers.flags['C_CONTIGUOUS']
        assert centers.shape[0] == ncharge
        self._charges = charges
        self._centers = centers
        self._this = <ints.GB2Integral*>(new ints.GB2GaussAttractionIntegral(
            max_nbasis, &charges[0], &centers[0, 0], ncharge, c, alpha
        ))

    property c:
        def __get__(self):
            return (<ints.GB2GaussAttractionIntegral*>self._this).get_c()

    property alpha:
        def __get__(self):
            return (<ints.GB2GaussAttractionIntegral*>self._this).get_alpha()


ints.libint2_static_init()
def libint2_static_cleanup():
    ints.libint2_static_cleanup()
atexit.register(libint2_static_cleanup)


cdef class GB4Integral:
    '''Wrapper for ints.GB4Integral'''
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
              np.ndarray[double, ndim=1] r0 not None, np.ndarray[double, ndim=1] r1 not None,
              np.ndarray[double, ndim=1] r2 not None, np.ndarray[double, ndim=1] r3 not None):
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
            np.ndarray[double, ndim=1] scales0 not None, np.ndarray[double, ndim=1] scales1 not None,
            np.ndarray[double, ndim=1] scales2 not None, np.ndarray[double, ndim=1] scales3 not None):
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


cdef class GB4ElectronRepulsionIntegralLibInt(GB4Integral):
    '''Wrapper for ints.GB4ElectronRepulsionIntegralLibInt, for testing only'''

    def __cinit__(self, long max_nbasis):
        self._this = <ints.GB4Integral*>(new ints.GB4ElectronRepulsionIntegralLibInt(max_nbasis))


cdef class GB4ErfIntegralLibInt(GB4Integral):
    '''Wrapper for ints.GB4ElectronRepulsionIntegralLibInt, for testing only'''

    def __cinit__(self, long max_nbasis, double mu):
        self._this = <ints.GB4Integral*>(new ints.GB4ErfIntegralLibInt(max_nbasis, mu))

    property mu:
        def __get__(self):
            return (<ints.GB4ErfIntegralLibInt*>self._this).get_mu()


cdef class GB4GaussIntegralLibInt(GB4Integral):
    '''Wrapper for ints.GB4GaussIntegralLibInt, for testing only'''

    def __cinit__(self, long max_nbasis, double c, double alpha):
        self._this = <ints.GB4Integral*>(new ints.GB4GaussIntegralLibInt(max_nbasis, c, alpha))

    property c:
        def __get__(self):
            return (<ints.GB4GaussIntegralLibInt*>self._this).get_c()

    property alpha:
        def __get__(self):
            return (<ints.GB4GaussIntegralLibInt*>self._this).get_alpha()


cdef class GB4RAlphaIntegralLibInt(GB4Integral):
    '''Wrapper for ints.GB4RAlphaIntegralLibInt, for testing only'''

    def __cinit__(self, long max_nbasis, double alpha):
        self._this = <ints.GB4Integral*>(new ints.GB4RAlphaIntegralLibInt(max_nbasis, alpha))

    property alpha:
        def __get__(self):
            return (<ints.GB4RAlphaIntegralLibInt*>self._this).get_alpha()


#
# fns wrappers (for testing and use in this module)
#


cdef class GB1DMGridFn:
    '''Wrapper for fns.GB1DMGridFn, for testing only'''
    cdef fns.GB1DMGridFn* _this

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

    property shell_type0:
        def __get__(self):
            return self._this.get_shell_type0()

    property dim_work:
        def __get__(self):
            return self._this.get_dim_work()

    property dim_output:
        def __get__(self):
            return self._this.get_dim_output()

    def reset(self, long shell_type0, np.ndarray[double, ndim=1] r0 not None, np.ndarray[double, ndim=1] point not None):
        assert r0.flags['C_CONTIGUOUS']
        assert r0.shape[0] == 3
        assert point.flags['C_CONTIGUOUS']
        assert point.shape[0] == 3
        self._this.reset(shell_type0, <double*>r0.data, &point[0])

    def add(self, double coeff, double alpha0,
            np.ndarray[double, ndim=1] scales0 not None):
        assert scales0.flags['C_CONTIGUOUS']
        assert scales0.shape[0] == get_shell_nbasis(abs(self._this.get_shell_type0()))
        self._this.add(coeff, alpha0, <double*>scales0.data)

    def cart_to_pure(self):
        self._this.cart_to_pure()

    def get_work(self, shape0):
        '''This returns a **copy** of the c++ work array.

           Returning a numpy array with a buffer created in c++ is dangerous.
           If the c++ array becomes deallocated, the numpy array may still
           point to the deallocated memory. For that reason, a copy is returned.
           Speed is not an issue as this class is only used for testing.
        '''
        cdef np.npy_intp shape[2]
        assert shape0 > 0
        assert shape0 <= self.max_nbasis
        shape[0] = shape0
        if self.dim_work == 1:
            tmp = np.PyArray_SimpleNewFromData(1, shape, np.NPY_DOUBLE, <void*> self._this.get_work())
        else:
            shape[1] = self.dim_work
            tmp = np.PyArray_SimpleNewFromData(2, shape, np.NPY_DOUBLE, <void*> self._this.get_work())
        return tmp.copy()


cdef class GB1DMGridDensityFn(GB1DMGridFn):
    def __cinit__(self, long max_nbasis):
        self._this = <fns.GB1DMGridFn*>(new fns.GB1DMGridDensityFn(max_nbasis))


cdef class GB1DMGridGradientFn(GB1DMGridFn):
    def __cinit__(self, long max_nbasis):
        self._this = <fns.GB1DMGridFn*>(new fns.GB1DMGridGradientFn(max_nbasis))


cdef class GB1DMGridGGAFn(GB1DMGridFn):
    def __cinit__(self, long max_nbasis):
        self._this = <fns.GB1DMGridFn*>(new fns.GB1DMGridGGAFn(max_nbasis))


cdef class GB1DMGridKineticFn(GB1DMGridFn):
    def __cinit__(self, long max_nbasis):
        self._this = <fns.GB1DMGridFn*>(new fns.GB1DMGridKineticFn(max_nbasis))


cdef class GB1DMGridHessianFn(GB1DMGridFn):
    def __cinit__(self, long max_nbasis):
        self._this = <fns.GB1DMGridFn*>(new fns.GB1DMGridHessianFn(max_nbasis))


cdef class GB1DMGridMGGAFn(GB1DMGridFn):
    def __cinit__(self, long max_nbasis):
        self._this = <fns.GB1DMGridFn*>(new fns.GB1DMGridMGGAFn(max_nbasis))


#
# iter_gb wrappers (for testing only)
#


cdef class IterGB1:
    """Wrapper for the IterGB1 class, for testing only."""
    cdef iter_gb.IterGB1* _this
    cdef GBasis _gbasis

    def __cinit__(self, GBasis gbasis not None):
        self._this = new iter_gb.IterGB1(gbasis._this)
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

    def store(self, np.ndarray[double, ndim=1] work not None,
              np.ndarray[double, ndim=1] output not None, long dim=1):
        max_shell_nbasis = get_shell_nbasis(self._gbasis.max_shell_type)
        assert work.shape[0] == get_shell_nbasis(self._this.shell_type0)
        assert work.flags['C_CONTIGUOUS']
        assert output.shape[0] == self._gbasis.nbasis
        assert output.flags['C_CONTIGUOUS']
        self._this.store(&work[0], &output[0], dim)

    property public_fields:
        def __get__(self):
            return (
                self._this.con_coeff,
                self._this.shell_type0,
                self._this.alpha0,
                self._this.r0[0], self._this.r0[1], self._this.r0[2],
                self._this.ibasis0,
            )

    property private_fields:
        def __get__(self):
            return (
                self._this.ishell0,
                self._this.nprim0,
                self._this.oprim0,
                self._this.iprim0,
            )


cdef class IterGB2:
    """Wrapper for the IterGB2 class, for testing only."""
    cdef iter_gb.IterGB2* _this
    cdef GBasis _gbasis

    def __cinit__(self, GBasis gbasis not None):
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

    def store(self, np.ndarray[double, ndim=2] work not None,
              np.ndarray[double, ndim=2] output not None):
        max_shell_nbasis = get_shell_nbasis(self._gbasis.max_shell_type)
        assert work.shape[0] == get_shell_nbasis(self._this.shell_type0)
        assert work.shape[1] == get_shell_nbasis(self._this.shell_type1)
        assert work.flags['C_CONTIGUOUS']
        assert output.shape[0] == self._gbasis.nbasis
        assert output.shape[1] == self._gbasis.nbasis
        assert output.flags['C_CONTIGUOUS']
        self._this.store(&work[0, 0], &output[0, 0])

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

    def __cinit__(self, GBasis gbasis not None):
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

    def store(self, np.ndarray[double, ndim=4] work not None,
              np.ndarray[double, ndim=4] output not None):
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
        self._this.store(&work[0, 0, 0, 0], &output[0, 0, 0, 0])

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


def iter_pow1_inc(np.ndarray[long, ndim=1] n not None):
    assert n.flags['C_CONTIGUOUS']
    assert n.shape[0] == 3
    return iter_pow.iter_pow1_inc(&n[0])


cdef class IterPow1:
    """Wrapper for the IterPow1 class, for testing only."""
    cdef iter_pow.IterPow1* _c_i1p

    def __cinit__(self):
        self._c_i1p = new iter_pow.IterPow1()

    def __dealloc__(self):
        del self._c_i1p

    def __init__(self, long shell_type0):
        if shell_type0 < 0:
            raise ValueError('A shell_type parameter can not be negative.')
        self._c_i1p.reset(shell_type0)

    def inc(self):
        return self._c_i1p.inc()

    property fields:
        def __get__(self):
            return (
                self._c_i1p.n0[0], self._c_i1p.n0[1], self._c_i1p.n0[2],
                self._c_i1p.ibasis0
            )


cdef class IterPow2:
    """Wrapper for the IterPow2 class, for testing only."""
    cdef iter_pow.IterPow2* _c_i2p

    def __cinit__(self):
        self._c_i2p = new iter_pow.IterPow2()

    def __dealloc__(self):
        del self._c_i2p

    def __init__(self, long shell_type0, long shell_type1):
        if shell_type0 < 0 or shell_type1 < 0:
            raise ValueError('A shell_type parameter can not be negative.')
        self._c_i2p.reset(shell_type0, shell_type1)

    def inc(self):
        return self._c_i2p.inc()

    property fields:
        def __get__(self):
            return (
                self._c_i2p.n0[0], self._c_i2p.n0[1], self._c_i2p.n0[2],
                self._c_i2p.n1[0], self._c_i2p.n1[1], self._c_i2p.n1[2],
                self._c_i2p.offset, self._c_i2p.ibasis0, self._c_i2p.ibasis1,
            )
