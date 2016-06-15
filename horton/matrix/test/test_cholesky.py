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


import numpy as np, h5py as h5
from nose.tools import assert_raises

from horton import *  # pylint: disable=wildcard-import,unused-wildcard-import


def test_linalg_factory_constructors():
    lf = CholeskyLinalgFactory(5)
    assert lf.default_nbasis == 5
    lf = CholeskyLinalgFactory()
    assert lf.default_nbasis is None
    lf.default_nbasis = 10

    # Four-index tests
    op4 = lf.create_four_index(nvec=8)
    assert isinstance(op4, CholeskyFourIndex)
    lf.create_four_index.__check_init_args__(lf, op4, nvec=8)
    assert op4.nbasis == 10
    assert op4.nvec == 8
    assert op4.shape == (10, 10, 10, 10)
    assert not op4.is_decoupled

    op4 = lf.create_four_index(8, 4)
    lf.create_four_index.__check_init_args__(lf, op4, 8, 4)
    assert op4.nbasis == 8
    assert op4.nvec == 4
    assert not op4.is_decoupled

    array = np.random.normal(0, 1, (5, 10, 10))
    op4 = lf.create_four_index(10, array=array)
    lf.create_four_index.__check_init_args__(lf, op4, nvec=5)
    assert op4._array is array
    assert op4._array2 is array
    assert op4.nbasis == 10
    assert op4.nvec == 5
    assert not op4.is_decoupled

    array2 = np.random.normal(0, 1, (5, 10, 10))
    op4 = lf.create_four_index(10, array=array, array2=array2)
    lf.create_four_index.__check_init_args__(lf, op4, nvec=5)
    assert op4._array is array
    assert op4._array2 is array2
    assert op4.nbasis == 10
    assert op4.nvec == 5
    assert op4.is_decoupled


def test_linalg_hdf5():
    # without default nbasis
    lf1 = CholeskyLinalgFactory()
    with h5.File('horton.matrix.test.test_cholesky.test_linalg_hdf5.h5', driver='core', backing_store=False) as f:
        lf1.to_hdf5(f)
        lf2 = CholeskyLinalgFactory.from_hdf5(f)
        assert isinstance(lf2, CholeskyLinalgFactory)
        assert lf2.default_nbasis is None
        lf3 = load_h5(f)
        assert isinstance(lf3, CholeskyLinalgFactory)
        assert lf3.default_nbasis is None

    # with default nbasis
    lf1 = CholeskyLinalgFactory(13)
    with h5.File('horton.matrix.test.test_cholesky.test_linalg_hdf5.h5', driver='core', backing_store=False) as f:
        lf1.to_hdf5(f)
        lf2 = CholeskyLinalgFactory.from_hdf5(f)
        assert isinstance(lf2, CholeskyLinalgFactory)
        assert lf2.default_nbasis == 13
        lf3 = load_h5(f)
        assert isinstance(lf3, CholeskyLinalgFactory)
        assert lf3.default_nbasis == 13


def test_linalg_objects_del():
    lf = CholeskyLinalgFactory()
    with assert_raises(TypeError):
        lf.create_four_index()


#
# Tests on the CholeskyFourIndex stuff
#


def get_four_cho_dense(nbasis=10, nvec=8, sym=8):
    '''Create random 2-index Cholesky vectors and matching dense four-index object

       **Optional arguments:**

       nbasis
            The number of basis functions

       nvec
            The number of Cholesky vectors

       sym
            The amount of symmetries in the FourIndex object. See
            :ref:`dense_matrix_symmetry` for more details.
    '''
    check_options('sym', sym, 1, 2, 4, 8)
    cho = CholeskyFourIndex(nbasis, nvec)
    if sym in (1, 4):
        cho.decouple_array2()
    cho.randomize()
    cho.symmetrize(sym)
    dense = DenseFourIndex(nbasis)
    dense._array[:] = np.einsum('kac,kbd->abcd', cho._array, cho._array2)
    assert dense.is_symmetric(sym)
    assert cho.is_symmetric(sym)
    return cho, dense


def test_four_index_hdf5():
    lf = CholeskyLinalgFactory(5)
    a = lf.create_four_index(5, 3)
    a.randomize()
    with h5.File('horton.matrix.test.test_cholesky.test_four_index_hdf5', driver='core', backing_store=False) as f:
        a.to_hdf5(f)
        b = CholeskyFourIndex.from_hdf5(f)
        assert a == b


def test_four_index_copy_new_randomize_clear_assign():
    lf = CholeskyLinalgFactory(5)
    for args in (None, 3), (4, 3):
        a = lf.create_four_index(*args)
        b = a.copy()
        b.randomize()
        assert a != b
        c = b.copy()
        c.new.__check_init_args__(c, b)
        assert b == c
        d = c.new()
        assert a == d
        b.assign(c)
        assert b == c


def test_four_index_iscale():
    lf = CholeskyLinalgFactory()
    op = lf.create_four_index(3, 2)
    op.randomize()
    tmp = op._array.copy()
    op.iscale(3.0)
    assert abs(op._array - 3**0.5*tmp).max() < 1e-10
    assert abs(op._array2 - 3**0.5*tmp).max() < 1e-10
    op.decouple_array2()
    op.randomize()
    tmp = op._array.copy()
    tmp2 = op._array2.copy()
    op.iscale(3.0)
    assert abs(op._array - 3**0.5*tmp).max() < 1e-10
    assert abs(op._array2 - 3**0.5*tmp2).max() < 1e-10


def test_four_index_get():
    cho, dense = get_four_cho_dense(nbasis=3, sym=1)
    for i0 in xrange(dense.shape[0]):
        for i1 in xrange(dense.shape[1]):
            for i2 in xrange(dense.shape[2]):
                for i3 in xrange(dense.shape[3]):
                    assert abs(cho.get_element(i0, i1, i2, i3) -
                               dense.get_element(i0, i1, i2, i3)) < 1e-10


def test_four_index_is_symmetric():
    for sym in 1, 2, 4, 8:
        cho = get_four_cho_dense(sym=sym)[0]
        assert cho.is_symmetric(sym)


def test_four_index_symmetrize():
    lf = CholeskyLinalgFactory(20)
    op = lf.create_four_index(nvec=8)
    for symmetry in 1, 2, 4, 8:
        op.decouple_array2()
        op.randomize()
        op.symmetrize(symmetry)
        assert op.is_symmetric(symmetry, 0, 0)


def test_four_index_symmetrize_order_of_operations():
    lf = CholeskyLinalgFactory(20)
    op = lf.create_four_index(nvec=8)
    for symmetry in 1, 2, 4, 8:
        op.decouple_array2()
        # ugly hack to have matrix elements with very different order of
        # magnitudes
        op._array[:] = 10**np.random.uniform(-20,20, (8,20,20))
        op._array2[:] = 10**np.random.uniform(-20,20, (8,20,20))
        op.symmetrize(symmetry)
        assert op.is_symmetric(symmetry, 0, 0)


def test_four_index_itranspose():
    for sym in 1, 2, 4, 8:
        cho = get_four_cho_dense(sym=sym)[0]
        cho.itranspose()
        assert cho.is_symmetric(sym)


def check_four_sum(sym):
    cho, dense = get_four_cho_dense(sym=sym)
    assert np.allclose(dense.sum(), cho.sum())


def test_four_sum_1():
    check_four_sum(1)


def test_four_sum_2():
    check_four_sum(2)


def test_four_sum_4():
    check_four_sum(4)


def test_four_sum_8():
    check_four_sum(8)


def check_four_slice_to_X(sym):
    cho, dense = get_four_cho_dense(sym=sym)

    for subscripts in "abab->ab", "aabb->ab", "abba->ab":
        # Return value
        factor = np.random.uniform(1, 2)
        assert np.allclose(dense.slice_to_two(subscripts, factor=factor)._array,
                           cho.slice_to_two(subscripts, factor=factor)._array)
        # Output argument
        dense_out = DenseTwoIndex(dense.nbasis)
        cho_out = DenseTwoIndex(cho.nbasis)
        dense.slice_to_two(subscripts, out=dense_out, factor=factor)
        cho.slice_to_two(subscripts, out=cho_out, factor=factor)
        assert np.allclose(dense_out._array, cho_out._array)
        # Output argument without clear
        factor = np.random.uniform(1, 2)
        dense.slice_to_two(subscripts, out=dense_out, factor=factor, clear=False)
        cho.slice_to_two(subscripts, out=cho_out, factor=factor, clear=False)
        assert np.allclose(dense_out._array, cho_out._array)

    for subscripts in 'abcc->bac', 'abcc->abc', 'abcb->abc', 'abbc->abc':
        # Return value
        factor = np.random.uniform(1, 2)
        assert np.allclose(dense.slice_to_three(subscripts, factor=factor)._array,
                           cho.slice_to_three(subscripts, factor=factor)._array)
        # Output argument
        dense_out = DenseThreeIndex(dense.nbasis)
        cho_out = DenseThreeIndex(cho.nbasis)
        dense.slice_to_three(subscripts, out=dense_out, factor=factor)
        cho.slice_to_three(subscripts, out=cho_out, factor=factor)
        assert np.allclose(dense_out._array, cho_out._array)
        # Output argument without clear
        factor = np.random.uniform(1, 2)
        dense.slice_to_three(subscripts, out=dense_out, factor=factor, clear=False)
        cho.slice_to_three(subscripts, out=cho_out, factor=factor, clear=False)
        assert np.allclose(dense_out._array, cho_out._array)


def test_four_slice_to_X_1():
    check_four_slice_to_X(1)


def test_four_slice_to_X_2():
    check_four_slice_to_X(2)


def test_four_slice_to_X_4():
    check_four_slice_to_X(4)


def test_four_slice_to_X_8():
    check_four_slice_to_X(8)


def check_four_contract_two_to_two_direct(sym):
    cho, dense = get_four_cho_dense(sym=sym)
    dm = DenseTwoIndex(dense.nbasis)
    dm.randomize()
    factor = np.random.uniform(1, 2)
    # check return value
    assert np.allclose(dense.contract_two_to_two('abcd,bd->ac', dm, factor=factor)._array,
                       cho.contract_two_to_two('abcd,bd->ac', dm, factor=factor)._array)
    # check output argument
    out_dense = DenseTwoIndex(dense.nbasis)
    out_cho = DenseTwoIndex(cho.nbasis)
    dense.contract_two_to_two('abcd,bd->ac', dm, out_dense, factor)
    cho.contract_two_to_two('abcd,bd->ac', dm, out_cho, factor)
    assert np.allclose(out_dense._array, out_cho._array)
    # check output argument abd clear=False
    factor = np.random.uniform(1, 2)
    dense.contract_two_to_two('abcd,bd->ac', dm, out_dense, factor, clear=False)
    cho.contract_two_to_two('abcd,bd->ac', dm, out_cho, factor, clear=False)
    assert np.allclose(out_dense._array, out_cho._array)



def test_four_contract_two_to_two_direct_1():
    check_four_contract_two_to_two_direct(1)


def test_four_contract_two_to_two_direct_2():
    check_four_contract_two_to_two_direct(2)


def test_four_contract_two_to_two_direct_4():
    check_four_contract_two_to_two_direct(4)


def test_four_contract_two_to_two_direct_8():
    check_four_contract_two_to_two_direct(8)


def check_four_contract_two_to_two_exchange(sym):
    cho, dense = get_four_cho_dense(sym=sym)
    dm = DenseTwoIndex(dense.nbasis)
    dm.randomize()
    factor = np.random.uniform(1, 2)
    # check return value
    assert np.allclose(dense.contract_two_to_two('abcd,cb->ad', dm, factor=factor)._array,
                       cho.contract_two_to_two('abcd,cb->ad', dm, factor=factor)._array)
    # check output argument
    out_dense = DenseTwoIndex(dense.nbasis)
    out_cho = DenseTwoIndex(cho.nbasis)
    dense.contract_two_to_two('abcd,cb->ad', dm, out_dense, factor)
    cho.contract_two_to_two('abcd,cb->ad', dm, out_cho, factor)
    assert np.allclose(out_dense._array, out_cho._array)
    # check output argument abd clear=False
    factor = np.random.uniform(1, 2)
    dense.contract_two_to_two('abcd,cb->ad', dm, out_dense, factor, clear=False)
    cho.contract_two_to_two('abcd,cb->ad', dm, out_cho, factor, clear=False)
    assert np.allclose(out_dense._array, out_cho._array)


def test_four_contract_two_to_two_exchange_1():
    check_four_contract_two_to_two_exchange(1)


def test_four_contract_two_to_two_exchange_2():
    check_four_contract_two_to_two_exchange(2)


def test_four_contract_two_to_two_exchange_4():
    check_four_contract_two_to_two_exchange(4)


def test_four_contract_two_to_two_exchange_8():
    check_four_contract_two_to_two_exchange(8)


def check_four_index_transform(sym_in, sym_exp, method):
    '''Test driver for four-index transform

       **Arguments:**

       sym_in
            The symmetry of the four-index object in the AO basis.

       sym_exp
            The symmetry of the orbitals used for the four-index transform.

       method
            'tensordot' or 'einsum'
    '''
    cho, dense = get_four_cho_dense(sym=sym_in)
    dense_mo = dense.new()
    cho_mo = cho.new()
    if sym_exp == 8:
        exp0 = DenseExpansion(dense.nbasis)
        exp0.randomize()
        dense_mo.assign_four_index_transform(dense, exp0, method=method)
        cho_mo.assign_four_index_transform(cho, exp0, method=method)
        assert cho_mo.is_decoupled == cho.is_decoupled
    elif sym_exp == 4:
        exp0 = DenseExpansion(dense.nbasis)
        exp0.randomize()
        exp1 = DenseExpansion(dense.nbasis)
        exp1.randomize()
        dense_mo.assign_four_index_transform(dense, exp0, exp1, method=method)
        cho_mo.assign_four_index_transform(cho, exp0, exp1, method=method)
        assert cho_mo.is_decoupled
    elif sym_exp == 2:
        exp0 = DenseExpansion(dense.nbasis)
        exp0.randomize()
        exp2 = DenseExpansion(dense.nbasis)
        exp2.randomize()
        dense_mo.assign_four_index_transform(dense, exp0, exp2=exp2, method=method)
        cho_mo.assign_four_index_transform(cho, exp0, exp2=exp2, method=method)
        assert cho_mo.is_decoupled == cho.is_decoupled
    elif sym_exp == 1:
        exp0 = DenseExpansion(dense.nbasis)
        exp0.randomize()
        exp1 = DenseExpansion(dense.nbasis)
        exp1.randomize()
        exp2 = DenseExpansion(dense.nbasis)
        exp2.randomize()
        exp3 = DenseExpansion(dense.nbasis)
        exp3.randomize()
        dense_mo.assign_four_index_transform(dense, exp0, exp1, exp2, exp3, method=method)
        cho_mo.assign_four_index_transform(cho, exp0, exp1, exp2, exp3, method=method)
        assert cho_mo.is_decoupled
    else:
        raise ValueError
    assert np.allclose(dense_mo._array, cho_mo.get_dense()._array)
    sym_and = symmetry_and(sym_in, sym_exp)
    assert cho_mo.is_symmetric(sym_and)
    assert dense_mo.is_symmetric(sym_and)


def symmetry_and(sym1, sym2):
    def to_mask(sym):
        return {1: 0, 2: 1, 4: 2, 8: 3}[sym]
    def from_mask(mask):
        return {0: 1, 1: 2, 2: 4, 3: 8}[mask]
    return from_mask(to_mask(sym1) & to_mask(sym2))


def test_four_index_transform_8_8_tensordot():
    check_four_index_transform(8, 8, 'tensordot')


def test_four_index_transform_8_8_einsum():
    check_four_index_transform(8, 8, 'einsum')


def test_four_index_transform_4_8_tensordot():
    check_four_index_transform(4, 8, 'tensordot')


def test_four_index_transform_4_8_einsum():
    check_four_index_transform(4, 8, 'einsum')


def test_four_index_transform_2_8_tensordot():
    check_four_index_transform(2, 8, 'tensordot')


def test_four_index_transform_2_8_einsum():
    check_four_index_transform(2, 8, 'einsum')


def test_four_index_transform_1_8_tensordot():
    check_four_index_transform(1, 8, 'tensordot')


def test_four_index_transform_1_8_einsum():
    check_four_index_transform(1, 8, 'einsum')


def test_four_index_transform_8_4_tensordot():
    check_four_index_transform(8, 4, 'tensordot')


def test_four_index_transform_8_4_einsum():
    check_four_index_transform(8, 4, 'einsum')


def test_four_index_transform_4_4_tensordot():
    check_four_index_transform(4, 4, 'tensordot')


def test_four_index_transform_4_4_einsum():
    check_four_index_transform(4, 4, 'einsum')


def test_four_index_transform_2_4_tensordot():
    check_four_index_transform(2, 4, 'tensordot')


def test_four_index_transform_2_4_einsum():
    check_four_index_transform(2, 4, 'einsum')


def test_four_index_transform_1_4_tensordot():
    check_four_index_transform(1, 4, 'tensordot')


def test_four_index_transform_1_4_einsum():
    check_four_index_transform(1, 4, 'einsum')


def test_four_index_transform_8_2_tensordot():
    check_four_index_transform(8, 2, 'tensordot')


def test_four_index_transform_8_2_einsum():
    check_four_index_transform(8, 2, 'einsum')


def test_four_index_transform_4_2_tensordot():
    check_four_index_transform(4, 2, 'tensordot')


def test_four_index_transform_4_2_einsum():
    check_four_index_transform(4, 2, 'einsum')


def test_four_index_transform_2_2_tensordot():
    check_four_index_transform(2, 2, 'tensordot')


def test_four_index_transform_2_2_einsum():
    check_four_index_transform(2, 2, 'einsum')


def test_four_index_transform_1_2_tensordot():
    check_four_index_transform(1, 2, 'tensordot')


def test_four_index_transform_1_2_einsum():
    check_four_index_transform(1, 2, 'einsum')


def test_four_index_transform_8_1_tensordot():
    check_four_index_transform(8, 1, 'tensordot')


def test_four_index_transform_8_1_einsum():
    check_four_index_transform(8, 1, 'einsum')


def test_four_index_transform_4_1_tensordot():
    check_four_index_transform(4, 1, 'tensordot')


def test_four_index_transform_4_1_einsum():
    check_four_index_transform(4, 1, 'einsum')


def test_four_index_transform_2_1_tensordot():
    check_four_index_transform(2, 1, 'tensordot')


def test_four_index_transform_2_1_einsum():
    check_four_index_transform(2, 1, 'einsum')


def test_four_index_transform_1_1_tensordot():
    check_four_index_transform(1, 1, 'tensordot')


def test_four_index_transform_1_1_einsum():
    check_four_index_transform(1, 1, 'einsum')
