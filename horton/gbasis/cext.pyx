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
#cython: embedsignature=True
'''C++ extensions'''


import numpy as np
cimport numpy as np
np.import_array()

cimport boys
cimport cartpure
cimport common
cimport gbasis
cimport ints
cimport fns
cimport iter_gb
cimport iter_pow

import atexit

from horton.log import log
from horton.matrix import LinalgFactory


__all__ = [
    # boys
    'boys_function',
    # cartpure
    'cart_to_pure_low',
    # common
    'fac', 'fac2', 'binom', 'get_shell_nbasis', 'get_max_shell_type',
    'gpt_coeff', 'gb_overlap_int1d', 'nuclear_attraction_helper',
    # gbasis
    'gob_cart_normalization', 'gob_pure_normalization',
    'GOBasis',
    # ints
    'GB2OverlapIntegral', 'GB2KineticIntegral',
    'GB2NuclearAttractionIntegral',
    'GB4ElectronReuplsionIntegralLibInt',
    # fns
    'GB1DMGridDensityFn', 'GB1DMGridGradientFn',
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
    def from_hdf5(cls, grp, lf):
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

    # low-level compute routines
    def compute_grid_point1(self, np.ndarray[double, ndim=1] output not None,
                            np.ndarray[double, ndim=1] point not None,
                            GB1DMGridFn grid_fn not None):
        assert output.flags['C_CONTIGUOUS']
        assert output.shape[0] == self.nbasis
        assert point.flags['C_CONTIGUOUS']
        assert point.shape[0] == 3
        self._this.compute_grid_point1(&output[0], &point[0], grid_fn._this)


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

    def compute_overlap(self, output):
        """Compute the overlap matrix in a Gaussian orbital basis.

           **Arguments:**

           output
                This can either be a OneBody instance (used to write the output
                to) or a LinalgFactory (used to allocate the output operator).
                In both cases, the resulting operator is returned.
        """
        # prepare the output array
        cdef np.ndarray[double, ndim=2] output_array
        if isinstance(output, LinalgFactory):
            lf = output
            output = lf.create_one_body(self.nbasis)
        output_array = output._array
        self.check_matrix_one_body(output_array)
        # call the low-level routine
        (<gbasis.GOBasis*>self._this).compute_overlap(&output_array[0, 0])
        # done
        return output

    def compute_kinetic(self, output):
        """Compute the kinetic energy matrix in a Gaussian orbital basis."""
        # prepare the output array
        cdef np.ndarray[double, ndim=2] output_array
        if isinstance(output, LinalgFactory):
            lf = output
            output = lf.create_one_body(self.nbasis)
        output_array = output._array
        self.check_matrix_one_body(output_array)
        # call the low-level routine
        (<gbasis.GOBasis*>self._this).compute_kinetic(&output_array[0, 0])
        # done
        return output

    def compute_nuclear_attraction(self,
                                   np.ndarray[double, ndim=1] charges not None,
                                   np.ndarray[double, ndim=2] centers not None,
                                   output):
        """Compute the kintic energy matrix in a Gaussian orbital basis."""
        # check arguments
        assert charges.flags['C_CONTIGUOUS']
        cdef long ncharge = charges.shape[0]
        assert centers.flags['C_CONTIGUOUS']
        assert centers.shape[0] == ncharge
        assert centers.shape[1] == 3
        # prepare the output array
        cdef np.ndarray[double, ndim=2] output_array
        if isinstance(output, LinalgFactory):
            lf = output
            output = lf.create_one_body(self.nbasis)
        output_array = output._array
        self.check_matrix_one_body(output_array)
        # call the low-level routine
        (<gbasis.GOBasis*>self._this).compute_nuclear_attraction(
            &charges[0], &centers[0, 0], ncharge,
            &output_array[0, 0],
        )
        # done
        return output

    def compute_electron_repulsion(self, output):
        # prepare the output array
        cdef np.ndarray[double, ndim=4] output_array
        if isinstance(output, LinalgFactory):
            lf = output
            output = lf.create_two_body(self.nbasis)
        output_array = output._array
        self.check_matrix_two_body(output_array)
        # call the low-level routine
        (<gbasis.GOBasis*>self._this).compute_electron_repulsion(&output_array[0, 0, 0, 0])
        # done
        return output

    def compute_grid_orbitals_exp(self, exp,
                                  np.ndarray[double, ndim=2] points not None,
                                  np.ndarray[long, ndim=1] iorbs not None,
                                  np.ndarray[double, ndim=2] orbs not None):
        '''Compute the orbtials on a grid for a given set of expansion coefficients.

           **Arguments:**

           exp
                An expansion object. For now, this must be a DenseExpansion object.

           points
                A Numpy array with grid points, shape (npoint,3).

           iorbs
                The indexes of the orbitals to be computed. If not given, the
                orbitals with a non-zero occupation number are computed

           orbs
                An output array, shape (npoint, len(iorbs)). The results are
                added to this array.

           **Warning:** the results are added to the output array!
        '''
        cdef np.ndarray[double, ndim=2] coeffs = exp.coeffs
        self.check_matrix_coeffs(coeffs)
        nfn = coeffs.shape[1]
        assert points.flags['C_CONTIGUOUS']
        npoint = points.shape[0]
        assert points.shape[1] == 3
        assert iorbs.flags['C_CONTIGUOUS']
        norb = iorbs.shape[0]
        assert orbs.flags['C_CONTIGUOUS']
        assert orbs.shape[0] == npoint
        assert orbs.shape[1] == norb
        (<gbasis.GOBasis*>self._this).compute_grid1_exp(
            nfn, &coeffs[0, 0], npoint, &points[0, 0],
            norb, &iorbs[0], &orbs[0, 0])

    def _compute_grid1_dm(self, dm, np.ndarray[double, ndim=2] points not None,
                          GB1DMGridFn grid_fn not None, np.ndarray output not None,
                          double epsilon=0):
        '''Compute some density function on a grid for a given density matrix.

           **Arguments:**

           dm
                A density matrix. For now, this must be a DenseOneBody object.

           points
                A Numpy array with grid points, shape (npoint,3).

           grid_fn
                A grid function.

           output
                A Numpy array for the output.

           **Optional arguments:**

           epsilon
                Allow errors on the density of this magnitude for the sake of
                efficiency.

           **Warning:** the results are added to the output array! This may
           be useful to combine results from different spin components.
        '''
        # Get the array of the density matrix
        cdef np.ndarray[double, ndim=2] dmar = dm._array
        self.check_matrix_one_body(dmar)

        # Get the maximum of the absolute value over the rows
        cdef np.ndarray[double, ndim=1] dmmaxrow = np.abs(dmar).max(axis=0)

        # Check the output array
        assert output.flags['C_CONTIGUOUS']
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
                                np.ndarray[double, ndim=1] rhos not None,
                                double epsilon=0):
        '''Compute the electron density on a grid for a given density matrix.

           **Arguments:**

           dm
                A density matrix. For now, this must be a DenseOneBody object.

           points
                A Numpy array with grid points, shape (npoint,3).

           rhos
                A Numpy array for the output, shape (npoint,).

           **Optional arguments:**

           epsilon
                Allow errors on the density of this magnitude for the sake of
                efficiency.

           **Warning:** the results are added to the output array! This may
           be useful to combine results from different spin components.
        '''
        self._compute_grid1_dm(dm, points, GB1DMGridDensityFn(self.max_shell_type), rhos, epsilon)

    def compute_grid_gradient_dm(self, dm,
                                 np.ndarray[double, ndim=2] points not None,
                                 np.ndarray[double, ndim=2] gradrhos not None,
                                 double epsilon=0):
        '''Compute the electron density gradient on a grid for a given density matrix.

           **Arguments:**

           dm
                A density matrix. For now, this must be a DenseOneBody object.

           points
                A Numpy array with grid points, shape (npoint,3).

           gradrhos
                A Numpy array for the output, shape (npoint,3).

           **Optional arguments:**

           epsilon
                Allow errors on the density of this magnitude for the sake of
                efficiency.

           **Warning:** the results are added to the output array! This may
           be useful to combine results from different spin components.
        '''
        self._compute_grid1_dm(dm, points, GB1DMGridGradientFn(self.max_shell_type), gradrhos, epsilon)

    def compute_grid_hartree_dm(self, dm,
                                np.ndarray[double, ndim=2] points not None,
                                np.ndarray[double, ndim=1] output not None):
        '''Compute the Hartree potential on a grid for a given density matrix.

           **Arguments:**

           dm
                A density matrix. For now, this must be a DenseOneBody object.

           points
                A Numpy array with grid points, shape (npoint,3).

           grid_fn
                A grid function.

           output
                A Numpy array for the output.

           **Warning:** the results are added to the output array! This may
           be useful to combine results from different spin components.
        '''
        cdef np.ndarray[double, ndim=2] dmar = dm._array
        self.check_matrix_one_body(dmar)
        assert output.flags['C_CONTIGUOUS']
        npoint = output.shape[0]
        assert points.flags['C_CONTIGUOUS']
        assert points.shape[0] == npoint
        assert points.shape[1] == 3
        (<gbasis.GOBasis*>self._this).compute_grid2_dm(
            &dmar[0, 0], npoint, &points[0, 0],
            &output[0])

    def _compute_grid1_fock(self, np.ndarray[double, ndim=2] points not None,
                           np.ndarray[double, ndim=1] weights not None,
                           np.ndarray pots not None,
                           GB1DMGridFn grid_fn not None, fock):
        '''Compute a one-body operator based on some potential grid in real-space

           **Arguments:**

           points
                A Numpy array with grid points, shape (npoint,3).

           weights
                A Numpy array with integration weights, shape (npoint,).

           pots
                A Numpy array with potential data on the grid.

           grid_fn
                A grid function.

           fock
                A one-body operator. For now, this must be a DenseOneBody
                object.

           **Warning:** the results are added to the fock operator!
        '''
        cdef np.ndarray[double, ndim=2] output = fock._array
        self.check_matrix_one_body(output)
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
        '''Compute a one-body operator based on a density potential grid in real-space

           **Arguments:**

           points
                A Numpy array with grid points, shape (npoint,3).

           weights
                A Numpy array with integration weights, shape (npoint,).

           pots
                A Numpy array with density potential data, shape (npoint,).

           fock
                A one-body operator. For now, this must be a DenseOneBody
                object.

           **Warning:** the results are added to the fock operator!
        '''
        self._compute_grid1_fock(points, weights, pots, GB1DMGridDensityFn(self.max_shell_type), fock)

    def compute_grid_gradient_fock(self, np.ndarray[double, ndim=2] points not None,
                                   np.ndarray[double, ndim=1] weights not None,
                                   np.ndarray[double, ndim=2] pots not None, fock):
        '''Compute a one-body operator based on a density potential grid in real-space

           **Arguments:**

           points
                A Numpy array with grid points, shape (npoint,3).

           weights
                A Numpy array with integration weights, shape (npoint,).

           pots
                A Numpy array with gradient potential data, shape (npoint, 3).

           fock
                A one-body operator. For now, this must be a DenseOneBody
                object.

           **Warning:** the results are added to the fock operator!
        '''
        self._compute_grid1_fock(points, weights, pots, GB1DMGridGradientFn(self.max_shell_type), fock)

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
    '''Wrapper for ints.GB2OverlapIntegral, for testing only'''

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


cdef class GB4ElectronReuplsionIntegralLibInt(GB4Integral):
    '''Wrapper for ints.GB4ElectronReuplsionIntegralLibInt, for testing only'''

    def __cinit__(self, long max_nbasis):
        self._this = <ints.GB4Integral*>(new ints.GB4ElectronReuplsionIntegralLibInt(max_nbasis))



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
