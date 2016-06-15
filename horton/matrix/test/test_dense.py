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


import numpy as np, h5py as h5, scipy as scipy
from nose.tools import assert_raises

from horton import *  # pylint: disable=wildcard-import,unused-wildcard-import


#
# Utility functions
#


def get_forth_back(n):
    '''Returns matching pair of forth and back permutation.

       **Arguments:**

       n
            The length of the permutation

       It is guaranteed that the identity permutation is never returned.
    '''
    while True:
        forth = np.random.uniform(0, 1, 5).argsort()
        if (forth != np.arange(5)).all():
            break
    back = np.zeros(5, int)
    for i, j in enumerate(forth):
        back[j] = i
    return forth, back


def get_signs(n):
    '''Returns an array with signs (all elements are just +1 or -1)

       **Arguments:**

       n
            The length of the signs array

       It is guaranteed that not all signs are positive.
    '''
    while True:
        signs = np.random.randint(0, 1, n)*2 -1
        if (signs < 0).any():
            return signs


def get_random_exp(lf):
    '''Return a random expansion and an identity overlap matrix'''
    exp = lf.create_expansion()
    a = np.random.normal(0, 1, (lf.default_nbasis, lf.default_nbasis))
    a = a + a.T
    evals, evecs = np.linalg.eigh(a)
    exp.coeffs[:] = evecs
    exp.occupations[:lf.default_nbasis/2] = 1.0
    exp.energies[:] = np.random.uniform(-1, 1, lf.default_nbasis)
    exp.energies.sort()
    olp = lf.create_two_index()
    olp._array[:] = np.identity(lf.default_nbasis)
    return exp, olp


#
# DenseLinalgFactory tests
#


def test_linalg_factory_constructors():
    lf = DenseLinalgFactory(5)
    assert lf.default_nbasis == 5
    lf = DenseLinalgFactory()
    assert lf.default_nbasis is None
    lf.default_nbasis = 10

    # One-index tests
    op1 = lf.create_one_index()
    assert isinstance(op1, DenseOneIndex)
    lf.create_one_index.__check_init_args__(lf, op1)
    assert op1.nbasis == 10
    assert op1.shape == (10,)
    op1 = lf.create_one_index(12)
    lf.create_one_index.__check_init_args__(lf, op1, 12)
    assert op1.nbasis == 12

    # Expansion tests
    ex = lf.create_expansion()
    assert isinstance(ex, DenseExpansion)
    lf.create_expansion.__check_init_args__(lf, ex)
    assert ex.nbasis == 10
    assert ex.nfn == 10
    assert ex.coeffs.shape == (10, 10)
    assert ex.energies.shape == (10,)
    assert ex.occupations.shape == (10,)
    ex = lf.create_expansion(12)
    lf.create_expansion.__check_init_args__(lf, ex, 12)
    assert ex.nbasis == 12
    assert ex.nfn == 12
    assert ex.coeffs.shape == (12, 12)
    assert ex.energies.shape == (12,)
    assert ex.occupations.shape == (12,)
    ex = lf.create_expansion(12, 10)
    lf.create_expansion.__check_init_args__(lf, ex, 12, 10)
    assert ex.nbasis == 12
    assert ex.nfn == 10
    assert ex.coeffs.shape == (12, 10)
    assert ex.energies.shape == (10,)
    assert ex.occupations.shape == (10,)

    # Two-index tests
    op2 = lf.create_two_index()
    assert isinstance(op2, DenseTwoIndex)
    lf.create_two_index.__check_init_args__(lf, op2)
    assert op2.nbasis == 10
    assert op2.shape == (10, 10)
    assert op2.nbasis1 == 10
    op2 = lf.create_two_index(12)
    lf.create_two_index.__check_init_args__(lf, op2, 12)
    assert op2.nbasis == 12
    assert op2.shape == (12, 12)
    assert op2.nbasis1 == 12
    op2 = lf.create_two_index(10, 12)
    lf.create_two_index.__check_init_args__(lf, op2, 10, 12)
    assert op2.shape == (10, 12)
    assert op2.nbasis == 10
    assert op2.nbasis1 == 12

    # Three-index tests
    op3 = lf.create_three_index()
    assert isinstance(op3, DenseThreeIndex)
    lf.create_three_index.__check_init_args__(lf, op3)
    assert op3.nbasis == 10
    assert op3.shape == (10, 10, 10)
    op3 = lf.create_three_index(8)
    lf.create_three_index.__check_init_args__(lf, op3, 8)
    assert op3.nbasis == 8

    # Four-index tests
    op4 = lf.create_four_index()
    assert isinstance(op4, DenseFourIndex)
    lf.create_four_index.__check_init_args__(lf, op4)
    assert op4.nbasis == 10
    assert op4.shape == (10, 10, 10, 10)
    op4 = lf.create_four_index(8)
    lf.create_four_index.__check_init_args__(lf, op4, 8)
    assert op4.nbasis == 8


def test_linalg_hdf5():
    # without default nbasis
    lf1 = DenseLinalgFactory()
    with h5.File('horton.matrix.test.test_dense.test_linalg_hdf5.h5', driver='core', backing_store=False) as f:
        lf1.to_hdf5(f)
        lf2 = DenseLinalgFactory.from_hdf5(f)
        assert isinstance(lf2, DenseLinalgFactory)
        assert lf2.default_nbasis is None
        lf3 = load_h5(f)
        assert isinstance(lf3, DenseLinalgFactory)
        assert lf3.default_nbasis is None

    # with default nbasis
    lf1 = DenseLinalgFactory(13)
    with h5.File('horton.matrix.test.test_dense.test_linalg_hdf5.h5', driver='core', backing_store=False) as f:
        lf1.to_hdf5(f)
        lf2 = DenseLinalgFactory.from_hdf5(f)
        assert isinstance(lf2, DenseLinalgFactory)
        assert lf2.default_nbasis == 13
        lf3 = load_h5(f)
        assert isinstance(lf3, DenseLinalgFactory)
        assert lf3.default_nbasis == 13


def test_linalg_objects_del():
    lf = DenseLinalgFactory()
    with assert_raises(TypeError):
        lf.create_one_index()
    with assert_raises(TypeError):
        lf.create_expansion()
    with assert_raises(TypeError):
        lf.create_two_index()
    with assert_raises(TypeError):
        lf.create_three_index()
    with assert_raises(TypeError):
        lf.create_four_index()


def test_allocate_check_output():
    # OneIndex
    lf = DenseLinalgFactory(5)
    op = lf.create_one_index()
    re = lf._allocate_check_output(op, (5,))
    assert op is re
    re = lf._allocate_check_output(None, (5,))
    assert isinstance(re, DenseOneIndex)
    assert re.shape == (5,)

    # TwoIndex
    lf = DenseLinalgFactory(5)
    op = lf.create_two_index()
    re = lf._allocate_check_output(op, (5, 5))
    assert op is re
    re = lf._allocate_check_output(None, (5, 5))
    assert isinstance(re, DenseTwoIndex)
    assert re.shape == (5, 5)

    # ThreeIndex
    lf = DenseLinalgFactory(5)
    op = lf.create_three_index()
    re = lf._allocate_check_output(op, (5, 5, 5))
    assert op is re
    re = lf._allocate_check_output(None, (5, 5, 5))
    assert isinstance(re, DenseThreeIndex)
    assert re.shape == (5, 5, 5)

    # FourIndex
    lf = DenseLinalgFactory(5)
    op = lf.create_four_index()
    re = lf._allocate_check_output(op, (5, 5, 5, 5))
    assert op is re
    re = lf._allocate_check_output(None, (5, 5, 5, 5))
    assert isinstance(re, DenseFourIndex)
    assert re.shape == (5, 5, 5, 5)


def test_einsum_wrapper():
    # all cases that occur in dense.py
    cases = [
        'ab->b', 'ab->a', 'ab,ab->b', 'ab,ab->a', 'ab,b->ab', 'ab,a->ab', 'aa',
        'ab', 'ab,ab', 'ab,ba', 'abc,ab->ac', 'abc,ab->ca', 'abc,bc->ba',
        'abc,bc->ab', 'abc,cb->ac', 'ab,c->cab', 'ab,c->acb', 'abc,db->adc',
        'abc,dc->adb', 'abc,db->dac', 'abc,dc->dab', 'aabb->ab', 'abab->ab',
        'abba->ab', 'abcc->bac', 'abcc->abc', 'abcb->abc', 'abbc->abc',
        'abcd,cd->acbd', 'abcd,cd->acdb', 'abcd,cb->acdb', 'abcd,cb->acbd',
        'abcd,ab->acbd', 'abcd,ab->acdb', 'abcd,ad->acbd', 'abcd,ad->acdb',
        'abcd,ad->abcd', 'abcd,ad->abdc', 'abcd,bd->abcd', 'abcd,bd->abdc',
        'abcd,bc->abdc', 'abcd,bc->abcd', 'abcd,ac->abcd', 'abcd,ac->abdc',
        'abcd,bd->ac', 'abcd,cb->ad'
    ]

    # run tests for every case
    lf = DenseLinalgFactory(5)
    for subscripts in cases:
        # prepare the inputs
        if '->' in subscripts:
            inscripts, outscript = subscripts.split('->')
        else:
            inscripts = subscripts
            outscript = None
        inscripts = inscripts.split(',')

        operands = [lf._allocate_check_output(None, (lf.default_nbasis,)*len(inscript)) for inscript in inscripts]
        for operand in operands:
            operand._array[:] = np.random.normal(0, 1, operand.shape)
        factor = np.random.uniform(1, 2)

        # compute with numpy directly
        if outscript is None:
            outarr = float(np.einsum(subscripts + '->...', *[operand._array for operand in operands]))*factor
        else:
            outarr = np.einsum(subscripts, *[operand._array for operand in operands])*factor

        # compute with wrapper and compare
        out = lf.einsum(subscripts, None, factor, True, *operands)

        if isinstance(out, float):
            assert abs(out - float(outarr)) < 1e-5
        else:
            assert outarr.shape == out.shape
            assert np.allclose(outarr, out._array)

        if not isinstance(out, float):
            # compute with wrapper (with output argument) and compare
            out = out.new()
            lf.einsum(subscripts, out, factor, True, *operands)
            assert np.allclose(outarr, out._array)

            # compute with wrapper (with output argument, without clear) and compare
            lf.einsum(subscripts, out, factor, False, *operands)
            assert np.allclose(outarr*2, out._array)


def test_tensordot_wrapper():
    # some cases to check
    cases = [
        (4, 2, ([1, 3], [1, 0])),
        (4, 2, ([1, 2], [1, 0])),
    ]

    lf = DenseLinalgFactory(5)
    for indim0, indim1, axes in cases:
        # construct inputs
        in0 = lf._allocate_check_output(None, (lf.default_nbasis,)*indim0)
        in0._array[:] = np.random.normal(0, 1, (in0.shape))
        in1 = lf._allocate_check_output(None, (lf.default_nbasis,)*indim1)
        in1._array[:] = np.random.normal(0, 1, (in1.shape))
        factor = np.random.uniform(1, 2)

        # get reference result as numpy array
        outarr = np.tensordot(in0._array, in1._array, axes)*factor

        # compute without output argument
        out = lf.tensordot(in0, in1, axes, None, factor)
        assert out.shape == outarr.shape
        assert np.allclose(out._array, outarr)

        # compute with output argument
        out = out.new()
        lf.tensordot(in0, in1, axes, out, factor)
        assert out.shape == outarr.shape
        assert np.allclose(out._array, outarr)

        # compute with output argument
        lf.tensordot(in0, in1, axes, out, factor, False)
        assert out.shape == outarr.shape
        assert np.allclose(out._array, outarr*2)


#
# DenseOneIndex tests
#


def test_one_index_hdf5():
    lf = DenseLinalgFactory(5)
    for args in (None,), (4,):
        a = lf.create_one_index(*args)
        a.randomize()
        with h5.File('horton.matrix.test.test_dens.test_one_index_hdf5', driver='core', backing_store=False) as f:
            a.to_hdf5(f)
            b = DenseOneIndex.from_hdf5(f)
            assert a == b


def test_one_index_copy_new_randomize_clear_assign():
    lf = DenseLinalgFactory(5)
    for args in (None,), (4,):
        a = lf.create_one_index(*args)
        b = a.copy()
        b.randomize()
        assert a != b
        c = b.copy()
        c.new.__check_init_args__(c, b)
        assert b == c
        d = c.new()
        assert a == d
        b.clear()
        assert a == b
        b.assign(c)
        assert b == c
        e = b.copy(1,3)
        assert e._array.shape[0] == 2
        assert (e._array == b._array[1:3]).all()
        b.clear()
        b.assign(e, 1, 3)
        assert (e._array == b._array[1:3]).all()
        assert ((e._array-b._array[1:3]) == 0).all()


def test_one_index_permute_basis():
    lf = DenseLinalgFactory(5)
    for i in xrange(10):
        forth, back = get_forth_back(5)
        a = lf.create_one_index()
        a.randomize()
        b = a.copy()
        b.permute_basis(forth)
        assert a != b
        b.permute_basis(back)
        assert a == b


def test_one_index_change_basis_signs():
    lf = DenseLinalgFactory(5)
    for i in xrange(10):
        signs = get_signs(5)
        a = lf.create_one_index()
        a.randomize()
        b = a.copy()
        b.change_basis_signs(signs)
        assert a != b
        b.change_basis_signs(signs)
        assert a == b


def test_one_index_iadd():
    lf = DenseLinalgFactory(5)
    a = lf.create_one_index()
    a.randomize()
    b = lf.create_one_index()
    b.randomize()
    c = b.copy()
    factor = np.random.uniform(1, 2)
    c.iadd(a, factor)
    for i in xrange(5):
        assert factor*a.get_element(i) + b.get_element(i) == c.get_element(i)


def test_one_index_iscale():
    lf = DenseLinalgFactory()
    op = lf.create_one_index(3)
    op.randomize()
    tmp = op._array.copy()
    op.iscale(3.0)
    assert abs(op._array - 3*tmp).max() < 1e-10


def test_one_index_norm():
    lf = DenseLinalgFactory(6)
    op = lf.create_one_index()
    op.randomize()
    norm = op.norm()
    assert (norm - np.linalg.norm(op._array)) < 1e-10


def test_one_index_get_set():
    lf = DenseLinalgFactory()
    op = lf.create_one_index(3)
    op.set_element(1, 1.2)
    assert op.get_element(1) == 1.2


def test_one_index_trace():
    lf = DenseLinalgFactory(5)
    inp = lf.create_one_index()
    inp.randomize()
    out = inp.trace()
    assert out == np.sum(inp._array)
    out = inp.trace(1, 3)
    assert out == np.sum(inp._array[1:3])


def test_one_index_get_max():
    lf = DenseLinalgFactory(5)
    inp = lf.create_one_index()
    inp.randomize()
    out = inp.get_max()
    assert out == np.max(np.abs(inp._array))


def test_one_index_sort_indices():
    lf = DenseLinalgFactory(5)
    op = lf.create_one_index()
    op.randomize()
    sortedlist = op.sort_indices(False)
    assert (sortedlist == np.argsort(op._array, axis=-1, kind='mergesort', order=None)[::-1]).all()
    sortedlistreverse = op.sort_indices(True)
    assert (sortedlistreverse == np.argsort(op._array, axis=-1, kind='mergesort', order=None)).all()


def test_one_index_mult():
    lf = DenseLinalgFactory(5)
    inp = lf.create_one_index()
    one = lf.create_one_index()
    inp.randomize()
    one.randomize()
    out = inp.mult(one, factor=1.3)
    assert np.allclose(out._array, 1.3*(inp._array*one._array))
    foo = inp.mult(one, out=out, factor=1.4)
    assert foo is out
    assert np.allclose(out._array, 1.4*(inp._array*one._array))


def test_one_index_dot():
    lf = DenseLinalgFactory(5)
    inp = lf.create_one_index()
    one = lf.create_one_index()
    inp.randomize()
    one.randomize()
    out = inp.dot(one, factor=1.3)
    assert out == 1.3*np.dot(inp._array, one._array)


def test_one_index_divide():
    lf = DenseLinalgFactory(5)
    inp = lf.create_one_index()
    one = lf.create_one_index()
    inp.randomize()
    one.randomize()
    out = inp.divide(one, factor=1.3)
    assert np.allclose(out._array, 1.3*np.divide(inp._array, one._array))
    foo = inp.divide(one, factor=1.4, out=out)
    assert foo is out
    assert np.allclose(out._array, 1.4*np.divide(inp._array, one._array))


#
# DenseExpansion tests
#


def test_expansion_hdf5():
    lf = DenseLinalgFactory(5)
    for args in (None,), (4,), (6, 3):
        a = lf.create_expansion(*args)
        a.randomize()
        with h5.File('horton.matrix.test.test_dens.test_expansion_hdf5', driver='core', backing_store=False) as f:
            a.to_hdf5(f)
            b = DenseExpansion.from_hdf5(f)
            assert a == b


def test_expansion_copy_new_randomize_clear_assign():
    lf = DenseLinalgFactory(5)
    for args in (None,), (4,), (6, 3):
        a = lf.create_expansion(*args)
        b = a.copy()
        b.randomize()
        assert a != b
        c = b.copy()
        c.new.__check_init_args__(c, b)
        assert b == c
        d = c.new()
        assert a == d
        b.clear()
        assert a == b
        b.assign(c)
        assert b == c


def test_expansion_copy():
    lf = DenseLinalgFactory()
    exp1 = lf.create_expansion(3, 2)
    exp1._coeffs[:] = np.random.uniform(0, 1, (3, 2))
    exp1._energies[:] = np.random.uniform(0, 1, 2)
    exp1._occupations[:] = np.random.uniform(0, 1, 2)
    exp2 = exp1.copy()
    assert (exp1._coeffs == exp2._coeffs).all()
    assert (exp1._energies == exp2._energies).all()
    assert (exp1._occupations == exp2._occupations).all()


def test_expansion_itranspose():
    lf = DenseLinalgFactory(5)
    orig = lf.create_expansion()
    one = lf.create_one_index()
    orig.randomize()
    out = orig.copy()
    out.itranspose()
    assert out.coeffs[0,1] == orig.coeffs[1,0]
    assert out.coeffs[2,1] == orig.coeffs[1,2]
    assert out.coeffs[3,2] == orig.coeffs[2,3]
    assert out.coeffs[4,2] == orig.coeffs[2,4]


def test_expansion_imul():
    lf = DenseLinalgFactory(5)
    orig = lf.create_expansion()
    one = lf.create_one_index()
    orig.randomize()
    one.randomize()
    out = orig.copy()
    out.imul(one)
    assert np.allclose(out.coeffs, orig.coeffs*one._array)


def test_expansion_permute_basis():
    lf = DenseLinalgFactory(5)
    for i in xrange(10):
        forth, back = get_forth_back(5)
        a = lf.create_expansion()
        a.randomize()
        b = a.copy()
        b.permute_basis(forth)
        assert a != b
        b.permute_basis(back)
        assert a == b


def test_expansion_permute_orbitals():
    lf = DenseLinalgFactory(5)
    for i in xrange(10):
        forth, back = get_forth_back(5)
        a = lf.create_expansion()
        a.randomize()
        b = a.copy()
        b.permute_orbitals(forth)
        assert a != b
        b.permute_orbitals(back)
        assert a == b


def test_expansion_change_basis_signs():
    lf = DenseLinalgFactory(5)
    for i in xrange(10):
        signs = get_signs(5)
        a = lf.create_expansion()
        a.randomize()
        b = a.copy()
        b.change_basis_signs(signs)
        assert a != b
        b.change_basis_signs(signs)
        assert a == b


def test_expansion_check_normalization():
    lf = DenseLinalgFactory(5)
    exp, olp = get_random_exp(lf)
    exp.check_normalization(olp)


def test_expansion_check_orthonormality():
    lf = DenseLinalgFactory(5)
    exp, olp = get_random_exp(lf)
    exp.check_orthonormality(olp)


def test_expansion_error_eigen():
    lf = DenseLinalgFactory(5)
    exp = lf.create_expansion()
    a = np.random.normal(0, 1, (5, 5))
    fock = lf.create_two_index()
    fock._array[:] = a+a.T
    evals, evecs = np.linalg.eigh(fock._array)
    exp.coeffs[:] = evecs
    exp.energies[:] = evals
    olp = lf.create_two_index()
    olp._array[:] = np.identity(5)
    assert exp.error_eigen(fock, olp) < 1e-10
    exp.coeffs[:] += np.random.normal(0, 1e-3, (5, 5))
    assert exp.error_eigen(fock, olp) > 1e-10


def test_expansion_from_fock():
    lf = DenseLinalgFactory(5)
    a = np.random.normal(0, 1, (5, 5))
    fock = lf.create_two_index()
    fock._array[:] = a+a.T
    a = np.random.normal(0, 1, (5, 5))
    olp = lf.create_two_index()
    olp._array[:] = np.dot(a, a.T)
    exp = lf.create_expansion()
    exp.from_fock(fock, olp)
    assert exp.error_eigen(fock, olp) < 1e-5


def test_expansion_from_fock_and_dm():
    natom = 5
    lf = DenseLinalgFactory(natom)

    # Use a simple Huckel-like model to construct degenerate levels
    fock = lf.create_two_index()
    olp = lf.create_two_index()
    for i in xrange(natom):
        fock.set_element(i, i, 0.6)
        fock.set_element(i, (i+1)%natom, -0.2)
        olp.set_element(i, i, 1.0)
        olp.set_element(i, (i+1)%natom, 0.2)

    # Create orbitals that will be used to construct various density matrices
    exp = lf.create_expansion()
    exp.from_fock(fock, olp)

    # Checks for every case
    def check_case(exp0):
        dm = lf.create_two_index()
        exp0.to_dm(dm)
        exp1 = lf.create_expansion()
        exp1.from_fock_and_dm(fock, dm, olp)
        assert np.allclose(exp0.occupations, exp1.occupations)
        assert exp1.error_eigen(fock, olp) < 1e-5
        sds = olp.copy()
        sds.itranspose()
        sds.idot(dm)
        sds.idot(olp)
        exp1.energies[:] = exp1.occupations
        assert exp1.error_eigen(sds, olp) < 1e-5

    # Case 1: not difficult, i.e. compatible degeneracies
    exp.occupations[:] = [1, 1, 1, 0, 0]
    check_case(exp)

    # Case 2: incompatible degeneracies
    exp.occupations[:] = [2, 2, 1, 0, 0]
    check_case(exp)

    # Case 3: incompatible degeneracies and rotated degenerate orbitals
    exp.occupations[:] = [2, 1, 0, 0, 0]
    for i in xrange(36):
        exp.rotate_2orbitals(np.pi/18.0, 1, 2)
        check_case(exp)

    # Case 4: incompatible degeneracies, fractional occupations and rotated
    # degenerate orbitals
    exp.occupations[:] = [1.5, 0.7, 0.3, 0, 0]
    for i in xrange(36):
        exp.rotate_2orbitals(np.pi/18.0, 1, 2)
        check_case(exp)


def test_expansion_naturals():
    fn_fchk = context.get_fn('test/ch3_hf_sto3g.fchk')
    mol = IOData.from_file(fn_fchk)
    overlap = mol.obasis.compute_overlap(mol.lf)
    dm = mol.exp_alpha.to_dm()
    exp = mol.lf.create_expansion()
    exp.derive_naturals(dm, overlap)
    assert exp.occupations.min() > -1e-6
    assert exp.occupations.max() < 1+1e-6
    exp.check_normalization(overlap)


def test_expansion_homo_lumo_ch3_hf():
    fn_fchk = context.get_fn('test/ch3_hf_sto3g.fchk')
    mol = IOData.from_file(fn_fchk)
    assert mol.exp_alpha.get_homo_index() == 4
    assert mol.exp_beta.get_homo_index() == 3
    assert mol.exp_alpha.get_lumo_index() == 5
    assert mol.exp_beta.get_lumo_index() == 4
    assert mol.exp_alpha.get_homo_index(1) == 3
    assert mol.exp_beta.get_homo_index(1) == 2
    assert mol.exp_alpha.get_lumo_index(1) == 6
    assert mol.exp_beta.get_lumo_index(1) == 5
    assert abs(mol.exp_alpha.get_homo_energy() - -3.63936540E-01) < 1e-8
    assert abs(mol.exp_alpha.get_homo_energy(1) - -5.37273275E-01) < 1e-8
    assert abs(mol.exp_alpha.get_lumo_energy() - 6.48361367E-01) < 1e-8
    assert abs(mol.exp_beta.get_homo_energy() - -5.18988806E-01) < 1e-8
    assert abs(mol.exp_beta.get_homo_energy(1) - -5.19454722E-01) < 1e-8
    assert abs(mol.exp_beta.get_lumo_energy() - 3.28562907E-01) < 1e-8
    assert abs(mol.exp_alpha.homo_energy - -3.63936540E-01) < 1e-8
    assert abs(mol.exp_alpha.lumo_energy - 6.48361367E-01) < 1e-8
    assert abs(mol.exp_beta.homo_energy - -5.18988806E-01) < 1e-8
    assert abs(mol.exp_beta.lumo_energy - 3.28562907E-01) < 1e-8


def test_expansion_to_dm1():
    fn_fchk = context.get_fn('test/water_hfs_321g.fchk')
    mol = IOData.from_file(fn_fchk)
    dm = mol.exp_alpha.to_dm()
    dm.iscale(2)
    assert dm.distance_inf(mol.get_dm_full()) < 1e-4
    assert dm.is_symmetric()


def test_expansion_to_dm2():
    fn_fchk = context.get_fn('test/ch3_hf_sto3g.fchk')
    mol = IOData.from_file(fn_fchk)
    dm = mol.exp_alpha.to_dm()
    dm1 = mol.exp_beta.to_dm(dm, 1.0, False)
    assert dm1 is dm
    olp = mol.obasis.compute_overlap(mol.lf)
    assert dm.distance_inf(mol.get_dm_full()) < 1e-4
    assert dm.is_symmetric()


def test_expansion_to_dm3():
    fn_fchk = context.get_fn('test/ch3_hf_sto3g.fchk')
    mol = IOData.from_file(fn_fchk)
    dm = mol.exp_alpha.to_dm(other=mol.exp_beta)
    assert not dm.is_symmetric()


def test_expansion_assign_dot1():
    lf = DenseLinalgFactory(2)
    exp0 = lf.create_expansion()
    exp1 = lf.create_expansion()
    tf2 = lf.create_two_index()
    exp0.randomize()
    tf2.randomize()
    exp1.assign_dot(exp0, tf2)
    assert np.allclose(exp1.coeffs, np.dot(exp0.coeffs, tf2._array))
    # exceptions
    exp3 = lf.create_expansion(nbasis=3, nfn=2)
    with assert_raises(TypeError):
        # mismatch between exp1.nbasis and exp3.nbasis
        exp1.assign_dot(exp3, tf2)
    with assert_raises(TypeError):
        # mismatch between exp3.nbasis and exp1.nbasis
        exp3.assign_dot(exp1, tf2)
    tf4 = lf.create_two_index(nbasis=3)
    exp5 = lf.create_expansion(nbasis=2, nfn=3)
    with assert_raises(TypeError):
        # mismatch between exp1.nfn and tf4.shape[0]
        exp5.assign_dot(exp1, tf4)
    exp4 = lf.create_expansion(nbasis=3, nfn=3)
    with assert_raises(TypeError):
        # mismatch between exp3.nfn and tf4.shape[1]
        exp3.assign_dot(exp4, tf4)


def test_expansion_assign_dot2():
    lf = DenseLinalgFactory(2)
    exp0 = lf.create_expansion()
    exp1 = lf.create_expansion()
    tf2 = lf.create_two_index()
    exp0.randomize()
    tf2.randomize()
    exp1.assign_dot(tf2, exp0)
    assert np.allclose(exp1.coeffs, np.dot(tf2._array, exp0.coeffs))
    # exceptions
    exp3 = lf.create_expansion(nbasis=3, nfn=2)
    with assert_raises(TypeError):
        # mismatch between tf2.shape[1] and exp3.nbasis
        exp1.assign_dot(tf2, exp3)
    with assert_raises(TypeError):
        # mismatch between tf2.shape[0] and exp3.nbasis
        exp3.assign_dot(tf2, exp1)
    tf4 = lf.create_two_index(nbasis=3)
    exp4 = lf.create_expansion(nbasis=3, nfn=3)
    exp5 = lf.create_expansion(nbasis=3, nfn=2)
    with assert_raises(TypeError):
        # mismatch between exp5.nfn and exp4.nfn
        exp5.assign_dot(tf4, exp4)


def test_expansion_assign_occupations():
    lf = DenseLinalgFactory(6)
    exp0 = lf.create_expansion()
    one = lf.create_one_index()
    exp0.randomize()
    one.randomize()
    exp0.assign_occupations(one)
    assert np.allclose(exp0.occupations, one._array)


def test_expansion_rotate_random():
    lf = DenseLinalgFactory(5)
    exp0, olp = get_random_exp(lf)
    exp0.check_normalization(olp)
    exp1 = exp0.copy()
    exp1.rotate_random()
    exp1.check_normalization(olp)
    dots = np.dot(exp0.coeffs.T, exp1.coeffs)
    assert not np.allclose(dots, np.identity(5))


def test_expansion_two_index_rotate_2orbitals():
    lf = DenseLinalgFactory(4)
    exp0, olp = get_random_exp(lf)
    exp0.check_normalization(olp)
    exp1 = exp0.copy()
    exp1.rotate_2orbitals()
    exp1.check_normalization(olp)
    check = np.identity(4, float)
    dots = np.dot(exp0.coeffs.T, exp1.coeffs)
    check = np.identity(4)
    check[1,1] = 1.0/np.sqrt(2)
    check[1,2] = 1.0/np.sqrt(2)
    check[2,1] = -1.0/np.sqrt(2)
    check[2,2] = 1.0/np.sqrt(2)
    assert np.allclose(dots, check)


def test_expansion_swap_orbitals():
    lf = DenseLinalgFactory(4)
    exp0, olp = get_random_exp(lf)
    exp0.check_normalization(olp)
    exp1 = exp0.copy()
    exp1.swap_orbitals(np.array([[0, 1], [2, 3]]))
    dots = np.dot(exp0.coeffs.T, exp1.coeffs)
    check = np.zeros((4, 4))
    check[0,1] = 1.0
    check[1,0] = 1.0
    check[2,3] = 1.0
    check[3,2] = 1.0
    assert np.allclose(dots, check)

#
# DenseTwoIndex tests
#


def test_two_index_hdf5():
    lf = DenseLinalgFactory(5)
    for args in (None,), (4,), (5, 6):
        a = lf.create_two_index(*args)
        a.randomize()
        with h5.File('horton.matrix.test.test_dens.test_two_index_hdf5', driver='core', backing_store=False) as f:
            a.to_hdf5(f)
            b = DenseTwoIndex.from_hdf5(f)
            assert a == b


def test_two_index_copy_new_randomize_clear_assign():
    lf = DenseLinalgFactory(5)
    for args in (None,), (4,), (5, 6):
        a = lf.create_two_index(*args)
        b = a.copy()
        b.randomize()
        assert a != b
        c = b.copy()
        c.new.__check_init_args__(c, b)
        assert b == c
        d = c.new()
        assert a == d
        b.clear()
        assert a == b
        b.assign(c)
        assert b == c
        e = b.copy(0, 3, 2, 4)
        assert e._array.shape[0] == 3
        assert e._array.shape[1] == 2
        assert (e._array == b._array[0:3, 2:4]).all()


def test_two_index_assign():
    lf = DenseLinalgFactory()
    op1 = lf.create_two_index(3)
    op2 = lf.create_two_index(3)
    op1._array[:] = np.random.uniform(0, 1, (3, 3))
    op2.assign(op1)
    assert (op1._array == op2._array).all()


def test_two_index_copy():
    lf = DenseLinalgFactory()
    op1 = lf.create_two_index(3)
    op1._array[:] = np.random.uniform(0, 1, (3, 3))
    op2 = op1.copy()
    assert (op1._array == op2._array).all()


def test_two_index_permute_basis():
    lf = DenseLinalgFactory(5)
    for i in xrange(10):
        i0, i1 = np.random.randint(0, 5, 2)
        forth, back = get_forth_back(5)
        a = lf.create_two_index()
        a.randomize()
        b = a.copy()
        b.permute_basis(forth)
        assert a != b
        assert a.get_element(i0, i1) == b.get_element(back[i0], back[i1])
        b.permute_basis(back)
        assert a == b


def test_two_index_change_basis_signs():
    lf = DenseLinalgFactory(5)
    for i in xrange(10):
        signs = get_signs(5)
        a = lf.create_two_index()
        a.randomize()
        b = a.copy()
        b.change_basis_signs(signs)
        assert a != b
        b.change_basis_signs(signs)
        assert a == b


def test_two_index_iadd():
    lf = DenseLinalgFactory()
    a = lf.create_two_index(5, 5)
    a.randomize()
    b = lf.create_two_index(5, 5)
    b.randomize()
    c = b.copy()
    factor = np.random.uniform(1, 2)
    # normal usage
    c.iadd(a, factor)
    for i0 in xrange(5):
        for i1 in xrange(5):
            assert factor*a.get_element(i0, i1) + b.get_element(i0, i1) == c.get_element(i0, i1)
    # transpose usage
    c.assign(b)
    c.iadd(a, factor, transpose=True)
    for i0 in xrange(5):
        for i1 in xrange(5):
            assert factor*a.get_element(i1, i0) + b.get_element(i0, i1) == c.get_element(i0, i1)
    # transpose usage
    c.assign(b)
    c.iadd_t(a, factor)
    for i0 in xrange(5):
        for i1 in xrange(5):
            assert factor*a.get_element(i1, i0) + b.get_element(i0, i1) == c.get_element(i0, i1)
    # constant
    c.assign(b)
    c.iadd(factor)
    for i0 in xrange(5):
        for i1 in xrange(5):
            assert factor + b.get_element(i0, i1) == c.get_element(i0, i1)
    # slice
    c.assign(b)
    a = lf.create_two_index(3, 3)
    a.randomize()
    c.iadd(a, factor, begin0=1, end0=4, begin1=1, end1=4)
    for i0 in xrange(5):
        for i1 in xrange(5):
            if i0 >= 1 and i0 < 4 and i1 >= 1 and i1 < 4:
                assert factor*a.get_element(i0-1, i1-1) + b.get_element(i0, i1) == c.get_element(i0, i1)
            else:
                assert b.get_element(i0, i1) == c.get_element(i0, i1)
    # slice and transpose
    c.assign(b)
    c.iadd(a, factor, begin0=1, end0=4, begin1=1, end1=4, transpose=True)
    for i0 in xrange(5):
        for i1 in xrange(5):
            if i0 >= 1 and i0 < 4 and i1 >= 1 and i1 < 4:
                assert factor*a.get_element(i1-1, i0-1) + b.get_element(i0, i1) == c.get_element(i0, i1)
            else:
                assert b.get_element(i0, i1) == c.get_element(i0, i1)


def test_two_index_iadd_slice():
    lf = DenseLinalgFactory()
    a = lf.create_two_index(5, 5)
    a.randomize()
    b = lf.create_two_index(8, 8)
    b.randomize()
    c = a.copy()
    factor = np.random.uniform(1, 2)
    # normal usage
    c.iadd_slice(b, factor, 1, 6, 3, 8)
    for i0 in xrange(5):
        for i1 in xrange(5):
            assert factor*b.get_element(i0+1, i1+3) + a.get_element(i0, i1) == c.get_element(i0, i1)


def test_two_index_iadd_one_mult():
    lf = DenseLinalgFactory()
    inp = lf.create_two_index(5, 5)
    inp.randomize()
    orig = inp.copy()
    one1 = lf.create_one_index(5)
    one1.randomize()
    one2 = lf.create_one_index(5)
    one2.randomize()
    factor = np.random.uniform(1, 2)
    # normal usage
    inp.iadd_one_mult(one1, one2, factor)
    for i0 in xrange(5):
        for i1 in xrange(5):
            assert abs(orig.get_element(i0, i1) + factor*one1.get_element(i1)*one2.get_element(i1) - inp.get_element(i0, i1)) < 1e-10
    # transpose usage
    inp.assign(orig)
    inp.iadd_one_mult(one1, one2, factor, transpose0=True)
    for i0 in xrange(5):
        for i1 in xrange(5):
            assert abs(orig.get_element(i0, i1) + factor*one1.get_element(i0)*one2.get_element(i1) - inp.get_element(i0, i1)) < 1e-10
    # transpose usage
    inp.assign(orig)
    inp.iadd_one_mult(one1, one2, factor, transpose1=True)
    for i0 in xrange(5):
        for i1 in xrange(5):
            assert abs(orig.get_element(i0, i1) + factor*one1.get_element(i1)*one2.get_element(i0) - inp.get_element(i0, i1)) < 1e-10
    # transpose usage
    inp.assign(orig)
    inp.iadd_one_mult(one1, one2, factor, transpose0=True, transpose1=True)
    for i0 in xrange(5):
        for i1 in xrange(5):
            assert abs(orig.get_element(i0, i1) + factor*one1.get_element(i0)*one2.get_element(i0) - inp.get_element(i0, i1)) < 1e-10


def test_two_index_iscale():
    lf = DenseLinalgFactory()
    op = lf.create_two_index(3)
    op.randomize()
    tmp = op._array.copy()
    op.iscale(3.0)
    assert abs(op._array - 3*tmp).max() < 1e-10


def test_two_index_get_set():
    lf = DenseLinalgFactory()
    op = lf.create_two_index(3)
    op.set_element(0, 1, 1.2)
    assert op.get_element(0, 1) == 1.2
    assert op.get_element(1, 0) == 1.2
    op = lf.create_two_index(3, 3)
    op.set_element(0, 1, 1.2, symmetry=1)
    assert op.get_element(0, 1) == 1.2
    assert op.get_element(1, 0) == 0.0


def test_two_index_sum():
    lf = DenseLinalgFactory()
    op1 = lf.create_two_index(3)
    op1._array[:] = np.random.uniform(-1, 1, (3,3))
    assert op1.sum() == op1._array.sum()
    assert op1.sum(begin0=1, end1=2) == op1._array[1,0] + op1._array[2,0] + op1._array[1,1] + op1._array[2,1]


def test_two_index_trace():
    lf = DenseLinalgFactory()
    op1 = lf.create_two_index(3)
    op1._array[:] = np.random.uniform(-1, 1, (3,3))
    assert op1.trace() == op1._array[0,0] + op1._array[1,1] + op1._array[2,2]
    assert op1.trace(begin0=1, end1=2) == op1._array[1,0] + op1._array[2,1]


def test_two_index_itranspose():
    lf = DenseLinalgFactory(8)
    for i in xrange(10):
        op = lf.create_two_index()
        op.randomize()
        i0, i1, = np.random.randint(0, 4, 2)
        x = op.get_element(i0, i1)
        op.itranspose()
        assert op.get_element(i1, i0) == x


def test_two_index_assign_dot():
    lf = DenseLinalgFactory(2)
    exp0 = lf.create_expansion()
    a = lf.create_two_index()
    tf2 = lf.create_two_index()
    exp0.randomize()
    tf2.randomize()
    a.assign_dot(exp0, tf2)
    assert np.allclose(a._array, np.dot(exp0.coeffs, tf2._array))


def test_two_index_inner():
    lf = DenseLinalgFactory(3)
    op = lf.create_two_index()
    vec0 = lf.create_one_index()
    vec1 = lf.create_one_index()
    op.randomize()
    vec0.randomize()
    vec1.randomize()
    assert abs(op.inner(vec0, vec1) - np.dot(vec0._array, np.dot(op._array, vec1._array))) < 1e-10
    assert abs(op.inner(vec0._array, vec1._array) - np.dot(vec0._array, np.dot(op._array, vec1._array))) < 1e-10


def test_two_index_sqrt():
    lf = DenseLinalgFactory(3)
    op = lf.create_two_index()
    op.randomize()
    op.imul(op)
    root = op.sqrt()
    assert np.allclose(root._array, scipy.linalg.sqrtm(op._array).real)


def test_two_index_inverse():
    a = np.array([[0.3,0.4,1.2],[3.0,1.2,0.5],[0.4,1.4,3.1]])
    lf = DenseLinalgFactory(3)
    op = lf.create_two_index()
    op.assign(a)
    unitm = op.new()
    unitm.assign_diagonal(1.0)
    inverse = op.inverse()
    assert np.allclose(np.dot(inverse._array, op._array), unitm._array)
    assert np.allclose(np.dot(op._array, inverse._array), unitm._array)


def test_two_index_assign_diagonal():
    lf = DenseLinalgFactory(3)
    op = lf.create_two_index()
    op.randomize()
    op.assign_diagonal(1.0)
    assert op.get_element(0, 0) == 1.0
    assert op.get_element(1, 1) == 1.0
    assert op.get_element(2, 2) == 1.0
    op.randomize()
    vec = lf.create_one_index()
    vec.randomize()
    op.assign_diagonal(vec)
    assert op.get_element(0, 0) == vec.get_element(0)
    assert op.get_element(1, 1) == vec.get_element(1)
    assert op.get_element(2, 2) == vec.get_element(2)


def test_two_index_copy_diagonal():
    lf = DenseLinalgFactory(3)
    op = lf.create_two_index()
    op.randomize()
    vec = op.copy_diagonal()
    assert op.get_element(0, 0) == vec.get_element(0)
    assert op.get_element(1, 1) == vec.get_element(1)
    assert op.get_element(2, 2) == vec.get_element(2)
    op.randomize()
    foo = op.copy_diagonal(vec)
    assert foo is vec
    assert op.get_element(0, 0) == vec.get_element(0)
    assert op.get_element(1, 1) == vec.get_element(1)
    assert op.get_element(2, 2) == vec.get_element(2)
    vec = op.copy_diagonal(begin=1)
    assert vec.shape == (2,)
    assert op.get_element(1, 1) == vec.get_element(0)
    assert op.get_element(2, 2) == vec.get_element(1)


def test_two_index_copy_slice():
    lf = DenseLinalgFactory(5)
    inp = lf.create_two_index()
    inp.randomize()
    ind = np.tril_indices(5, -1)
    out = inp.copy_slice(ind)
    assert inp.get_element(1, 0) == out.get_element(0)
    assert inp.get_element(2, 0) == out.get_element(1)
    assert inp.get_element(2, 1) == out.get_element(2)
    assert inp.get_element(3, 0) == out.get_element(3)
    assert inp.get_element(3, 1) == out.get_element(4)
    assert inp.get_element(3, 2) == out.get_element(5)
    assert inp.get_element(4, 0) == out.get_element(6)
    assert inp.get_element(4, 1) == out.get_element(7)
    assert inp.get_element(4, 2) == out.get_element(8)
    assert inp.get_element(4, 3) == out.get_element(9)
    foo = inp.copy_diagonal()
    foo = inp.copy_slice(ind, out)
    assert foo is out
    assert foo.get_element(0) == out.get_element(0)
    assert foo.get_element(1) == out.get_element(1)
    assert foo.get_element(2) == out.get_element(2)
    assert foo.get_element(3) == out.get_element(3)
    assert foo.get_element(4) == out.get_element(4)
    assert foo.get_element(5) == out.get_element(5)
    assert foo.get_element(6) == out.get_element(6)
    assert foo.get_element(7) == out.get_element(7)
    assert foo.get_element(8) == out.get_element(8)
    assert foo.get_element(9) == out.get_element(9)


def test_two_index_is_symmetric():
    lf = DenseLinalgFactory(4)
    op = lf.create_two_index()
    op.set_element(2, 3, 3.1234)
    assert op.is_symmetric()
    op.set_element(2, 3, 3.1, symmetry=1)
    assert not op.is_symmetric()
    assert op.is_symmetric(1)


def test_two_index_symmetrize():
    lf = DenseLinalgFactory(3)
    op = lf.create_two_index()
    op.randomize()
    assert not op.is_symmetric()
    op.symmetrize()
    x = op.get_element(1,2)
    op.set_element(2, 0, 0.0, symmetry=1)
    op.set_element(0, 2, 1.0, symmetry=1)
    op.symmetrize()
    op.is_symmetric()
    assert op.get_element(1,2) == x
    assert op.get_element(0,2) == 0.5


def test_two_index_contract_to_one():
    lf = DenseLinalgFactory(5)
    op = lf.create_two_index()
    op.randomize()
    op.symmetrize()
    # regular use
    vec = op.contract_to_one('ab->a')
    assert np.allclose(vec._array, op._array.sum(axis=0))
    vec = op.contract_to_one('ab->b')
    assert np.allclose(vec._array, op._array.sum(axis=1))
    # with ranges
    vec = op.contract_to_one('ab->b', end1=4, begin0=1)
    assert vec.shape == (4,)
    assert np.allclose(vec._array, op._array[:4,1:].sum(axis=1))
    vec = op.contract_to_one('ab->a', begin1=1, end0=4)
    assert vec.shape == (4,)
    assert np.allclose(vec._array, op._array[1:,:4].sum(axis=0))
    # with output argument, clear=True
    factor = np.random.uniform(1, 2)
    vec = lf.create_one_index()
    foo = op.contract_to_one('ab->a', vec, factor)
    assert foo is vec
    assert np.allclose(vec._array, factor*op._array.sum(axis=0))
    op.contract_to_one('ab->b', vec, factor)
    assert np.allclose(vec._array, factor*op._array.sum(axis=1))
    # with output argument, clear=False
    factor = np.random.uniform(1, 2)
    orig = lf.create_one_index()
    orig.randomize()
    vec = orig.copy()
    foo = op.contract_to_one('ab->a', vec, factor, clear=False)
    assert np.allclose(vec._array, orig._array + factor*op._array.sum(axis=0))
    vec = orig.copy()
    op.contract_to_one('ab->b', vec, factor, clear=False)
    assert np.allclose(vec._array, orig._array + factor*op._array.sum(axis=1))


def test_two_index_contract_two_to_one():
    lf = DenseLinalgFactory(5)
    a = lf.create_two_index()
    b = lf.create_two_index()
    c = lf.create_two_index(4, 4)
    a.randomize()
    a.symmetrize()
    b.randomize()
    b.symmetrize()
    c.randomize()
    c.symmetrize()
    # regular use
    vec = a.contract_two_to_one('ab,ab->a', b)
    assert np.allclose(vec._array, (a._array*b._array).sum(axis=1))
    vec = a.contract_two_to_one('ab,ab->b', b)
    assert np.allclose(vec._array, (a._array*b._array).sum(axis=0))
    # with output, clear=True
    factor = np.random.uniform(1, 2)
    foo = a.contract_two_to_one('ab,ab->a', b, vec, factor)
    assert foo is vec
    assert np.allclose(vec._array, factor*(a._array*b._array).sum(axis=1))
    a.contract_two_to_one('ab,ab->b', b, vec, factor)
    assert np.allclose(vec._array, factor*(a._array*b._array).sum(axis=0))
    # with output, clear=False
    factor = np.random.uniform(1, 2)
    orig = lf.create_one_index()
    orig.randomize()
    vec = orig.copy()
    a.contract_two_to_one('ab,ab->a', b, vec, factor, False)
    assert np.allclose(vec._array, orig._array + factor*(a._array*b._array).sum(axis=1))
    vec = orig.copy()
    a.contract_two_to_one('ab,ab->b', b, vec, factor, False)
    assert np.allclose(vec._array, orig._array + factor*(a._array*b._array).sum(axis=0))
    # with ranges
    vec = c.contract_two_to_one('ab,ab->a', b, end1=4, begin0=1)
    assert vec.shape == (4,)
    assert np.allclose(vec._array, (c._array*b._array[1:,:4]).sum(axis=1))
    vec = c.contract_two_to_one('ab,ab->b', b, end1=4, begin0=1)
    assert vec.shape == (4,)
    assert np.allclose(vec._array, (c._array*b._array[1:,:4]).sum(axis=0))
    # with ranges, output, clear=True
    factor = np.random.uniform(1, 2)
    foo = c.contract_two_to_one('ab,ab->a', b, vec, factor, end1=4, begin0=1)
    assert foo is vec
    assert np.allclose(vec._array, factor*(c._array*b._array[1:,:4]).sum(axis=1))
    c.contract_two_to_one('ab,ab->b', b, vec, factor, end1=4, begin0=1)
    assert np.allclose(vec._array, factor*(c._array*b._array[1:,:4]).sum(axis=0))
    # with ranges, output, clear=False
    factor = np.random.uniform(1, 2)
    orig = lf.create_one_index(4)
    orig.randomize()
    vec = orig.copy()
    c.contract_two_to_one('ab,ab->a', b, vec, factor, False, end1=4, begin0=1)
    assert np.allclose(vec._array, orig._array + factor*(c._array*b._array[1:,:4]).sum(axis=1))
    vec = orig.copy()
    c.contract_two_to_one('ab,ab->b', b, vec, factor, False, end1=4, begin0=1)
    assert np.allclose(vec._array, orig._array + factor*(c._array*b._array[1:,:4]).sum(axis=0))


def test_two_index_iadd_outer():
    lf = DenseLinalgFactory(5)
    a = lf.create_two_index()
    b = lf.create_two_index()
    a.randomize()
    b.randomize()
    orig = lf.create_two_index(5*5)
    orig.randomize()
    out = orig.copy()
    out.iadd_outer(a, b, 2.5)
    assert abs(out.get_element(0, 0) - orig.get_element(0, 0) - 2.5*a.get_element(0, 0)*b.get_element(0, 0)) < 1e-8
    assert abs(out.get_element(18, 11) - orig.get_element(18, 11) - 2.5*a.get_element(3, 3)*b.get_element(2, 1)) < 1e-8


def test_two_index_iadd_kron():
    lf = DenseLinalgFactory(5)
    a = lf.create_two_index()
    b = lf.create_two_index()
    a.randomize()
    b.randomize()
    orig = lf.create_two_index(5*5)
    orig.randomize()
    out = orig.copy()
    out.iadd_kron(a, b, 2.5)
    assert abs(out.get_element(0, 0) - orig.get_element(0, 0) - 2.5*a.get_element(0, 0)*b.get_element(0, 0)) < 1e-8
    assert abs(out.get_element(18, 11) - orig.get_element(18, 11) - 2.5*a.get_element(3, 2)*b.get_element(3, 1)) < 1e-8


def test_two_index_iadd_dot_all():
    lf = DenseLinalgFactory(5)
    a = lf.create_two_index()
    b = lf.create_two_index()
    a.randomize()
    b.randomize()
    orig = lf.create_two_index()
    orig.randomize()
    out = orig.copy()
    out.iadd_dot(a, b)
    assert np.allclose(out._array, orig._array + np.dot(a._array, b._array))
    out = orig.copy()
    out.iadd_tdot(a, b)
    assert np.allclose(out._array, orig._array + np.dot(a._array.T, b._array))
    out = orig.copy()
    out.iadd_dott(a, b)
    assert np.allclose(out._array, orig._array + np.dot(a._array, b._array.T))
    # with ranges
    a = lf.create_two_index(7)
    b = lf.create_two_index(9)
    a.randomize()
    b.randomize()
    out = orig.copy()
    out.iadd_dot(a, b, begin1=2, end1=7, begin0=0, end0=5, begin2=3, end2=8, begin3=4, end3=9)
    assert np.allclose(out._array, orig._array + np.dot(a._array[:5,2:7], b._array[3:8,4:9]))


def test_two_index_iadd_mult():
    lf = DenseLinalgFactory(5)
    a = lf.create_two_index()
    b = lf.create_two_index()
    a.randomize()
    b.randomize()
    orig = lf.create_two_index()
    orig.randomize()
    out = orig.copy()
    out.iadd_mult(a, b)
    assert np.allclose(out._array, orig._array + (a._array*b._array))
    out = orig.copy()
    out.iadd_mult(a, b, transpose0=True)
    assert np.allclose(out._array, orig._array + (a._array.T*b._array))
    out = orig.copy()
    out.iadd_mult(a, b, transpose1=True)
    assert np.allclose(out._array, orig._array + (a._array*b._array.T))
    out = orig.copy()
    out.iadd_mult(a, b, transpose0=True, transpose1=True)
    assert np.allclose(out._array, orig._array + (a._array.T*b._array.T))
    # with ranges
    a = lf.create_two_index(7)
    b = lf.create_two_index(9)
    a.randomize()
    b.randomize()
    out = orig.copy()
    out.iadd_mult(b, a, begin1=2, end1=7, begin0=0, end0=5, begin2=3, end2=8, begin3=4, end3=9)
    assert np.allclose(out._array, orig._array + (b._array[3:8,4:9]*a._array[:5,2:7]))
    out = orig.copy()
    out.iadd_mult(b, a, begin1=2, end1=7, begin0=0, end0=5, begin2=3, end2=8, begin3=4, end3=9, transpose0=True)
    assert np.allclose(out._array, orig._array + (b._array[3:8,4:9].T*a._array[:5,2:7]))
    out = orig.copy()
    out.iadd_mult(b, a, begin1=2, end1=7, begin0=0, end0=5, begin2=3, end2=8, begin3=4, end3=9, transpose1=True)
    assert np.allclose(out._array, orig._array + (b._array[3:8,4:9]*a._array[:5,2:7].T))
    out = orig.copy()
    out.iadd_mult(b, a, begin1=2, end1=7, begin0=0, end0=5, begin2=3, end2=8, begin3=4, end3=9, transpose1=True, transpose0=True)
    assert np.allclose(out._array, orig._array + (b._array[3:8,4:9].T*a._array[:5,2:7].T))


def test_two_index_iadd_shift():
    lf = DenseLinalgFactory(5)
    orig = lf.create_two_index()
    orig.randomize()
    op = orig.copy()
    shift = 0.3
    op.iadd_shift(shift)
    assert abs(op._array.min()) >= shift
    for i in xrange(5):
        for j in xrange(5):
            if orig.get_element(i, j) < 0.0:
                assert op.get_element(i, j) == orig.get_element(i, j) - shift
            else:
                assert op.get_element(i, j) == orig.get_element(i, j) + shift


def test_two_index_iadd_contract_two_one():
    lf = DenseLinalgFactory(5)
    x = lf.create_two_index()
    y = lf.create_one_index()
    orig = lf.create_two_index()
    x.randomize()
    y.randomize()
    orig.randomize()
    z = orig.copy()
    z.iadd_contract_two_one('ab,b->ab', x, y, 3.4)
    assert np.allclose(z._array, orig._array + 3.4*x._array*y._array)
    z = orig.copy()
    z.iadd_contract_two_one('ab,a->ab', x, y, 3.4)
    assert np.allclose(z._array, orig._array + 3.4*x._array*y._array[:,None])
    z = orig.copy()
    z.iadd_contract_two_one('ba,a->ab', x, y, 3.4)
    assert np.allclose(z._array, orig._array + 3.4*x._array.T*y._array[:,None])
    # with ranges
    x = lf.create_two_index(7)
    y = lf.create_one_index(9)
    x.randomize()
    y.randomize()
    z = orig.copy()
    z.iadd_contract_two_one('ab,b->ab', x, y, 3.4, 2, 7, 0, 5, 3, 8)
    assert np.allclose(z._array, orig._array + 3.4*x._array[2:7,:5]*y._array[3:8])
    z = orig.copy()
    z.iadd_contract_two_one('ab,a->ab', x, y, 3.4, 2, 7, 0, 5, 3, 8)
    assert np.allclose(z._array, orig._array + 3.4*x._array[2:7,:5]*y._array[3:8,None])
    z = orig.copy()
    z.iadd_contract_two_one('ba,a->ab', x, y, 3.4, 2, 7, 0, 5, 3, 8)
    assert np.allclose(z._array, orig._array + 3.4*x._array[2:7,:5].T*y._array[3:8,None])


def test_two_index_contract_two_to_two():
    # Test in detail 'ab,cb->ac'
    lf = DenseLinalgFactory(5)
    inp = lf.create_two_index()
    two = lf.create_two_index()
    inp.randomize()
    two.randomize()
    out = inp.contract_two_to_two('ab,cb->ac', two, factor=1.3)
    assert np.allclose(out._array, 1.3*np.einsum('ab,cb->ac', inp._array, two._array))
    foo = inp.contract_two_to_two('ab,cb->ac', two, out, factor=1.4)
    assert foo is out
    assert np.allclose(out._array, 1.4*np.einsum('ab,cb->ac', inp._array, two._array))
    inp.contract_two_to_two('ab,cb->ac', two, out, factor=1.4, clear=False)
    assert np.allclose(out._array, 2.8*np.einsum('ab,cb->ac', inp._array, two._array))
    # Blind test on all other cases
    others = ['ab,ca->cb', 'ab,bc->ac', 'ab,ac->cb', 'ab,ac->bc', 'ab,cb->ca']
    for subscripts in others:
        inp.contract_two_to_two(subscripts, two, factor=1.7)
    # with ranges
    two = lf.create_two_index(9)
    two.randomize()
    out = inp.contract_two_to_two('ab,cb->ac', two, factor=1.3, begin1=1, end1=6, begin0=3, end0=8)
    assert np.allclose(out._array, 1.3*np.einsum('ab,cb->ac', inp._array, two._array[3:8, 1:6]))
    foo = inp.contract_two_to_two('ab,cb->ac', two, out, factor=1.4, begin1=1, end1=6, begin0=3, end0=8)
    assert foo is out
    assert np.allclose(out._array, 1.4*np.einsum('ab,cb->ac', inp._array, two._array[3:8, 1:6]))
    # Blind test on all other cases
    others = ['ab,ca->cb', 'ab,bc->ac', 'ab,ac->cb', 'ab,ac->bc', 'ab,cb->ca']
    for subscripts in others:
        inp.contract_two_to_two(subscripts, two, factor=1.7, begin1=1, end1=6, begin0=3, end0=8)


def test_two_index_contract_two_abba():
    lf = DenseLinalgFactory()
    op1 = lf.create_two_index(3)
    op2 = lf.create_two_index(3)
    op1._array[:] = np.random.uniform(0, 1, (3, 3))
    op2._array[:] = np.random.uniform(0, 1, (3, 3))

    value = op1.contract_two('ab,ba', op2)
    op1.idot(op2)
    assert abs(op1.trace() - value) < 1e-10
    # with ranges
    op2 = lf.create_two_index(5)
    op2._array[:] = np.random.uniform(0, 1, (5, 5))
    value = op1.contract_two('ab,ba', op2, 0, 3, 2, 5)
    assert abs(np.dot(op1._array, op2._array[:3, 2:5]).trace() - value) < 1e-10


def test_two_index_contract_two_abab():
    lf = DenseLinalgFactory()
    op1 = lf.create_two_index(3)
    op2 = lf.create_two_index(3)
    op1._array[:] = np.random.uniform(0, 1, (3, 3))
    op2._array[:] = np.random.uniform(0, 1, (3, 3))

    value = op1.contract_two('ab,ab', op2)
    op2.itranspose()
    op1.idot(op2)
    assert abs(op1.trace() - value) < 1e-10
    # with ranges
    op2 = lf.create_two_index(5)
    op2._array[:] = np.random.uniform(0, 1, (5, 5))
    value = op1.contract_two('ab,ab', op2, 0, 3, 2, 5)
    assert abs(np.dot(op1._array.ravel(), op2._array[:3, 2:5].ravel()) - value) < 1e-10


def test_two_index_idot():
    lf = DenseLinalgFactory(3)
    op1 = lf.create_two_index()
    op2 = lf.create_two_index()
    op1.randomize()
    op2.randomize()
    op3 = op2.copy()
    op2.idot(op1)
    assert np.allclose(op2._array, np.dot(op3._array, op1._array))


def test_two_index_imul():
    lf = DenseLinalgFactory(3)
    orig = lf.create_two_index()
    op2 = lf.create_two_index()
    orig.randomize()
    op2.randomize()
    op1 = orig.copy()
    op1.imul(op2, 1.3)
    assert np.allclose(op1._array, orig._array*op2._array*1.3)
    op1 = orig.copy()
    op1.imul_t(op2, 1.3)
    assert np.allclose(op1._array, orig._array*op2._array.T*1.3)
    # with ranges
    op2 = lf.create_two_index(5)
    op2.randomize()
    op1 = orig.copy()
    op1.imul(op2, 1.3, 0, 3, 2, 5)
    assert np.allclose(op1._array, orig._array*op2._array[:3,2:5]*1.3)
    op1 = orig.copy()
    op1.imul_t(op2, 1.3, 0, 3, 2, 5)
    assert np.allclose(op1._array, orig._array*op2._array[:3,2:5].T*1.3)


def test_two_index_distance_inf():
    lf = DenseLinalgFactory(3)
    a = lf.create_two_index()
    a.randomize()
    b = a.copy()
    assert np.isclose(a.distance_inf(b), 0.0)
    b.set_element(0, 0, b.get_element(0, 0)+0.1)
    assert np.isclose(a.distance_inf(b), 0.1)


def test_two_index_iabs():
    lf = DenseLinalgFactory(3)
    a = lf.create_two_index()
    a.randomize()
    b = a.copy()
    a.iabs()
    assert np.allclose(a._array, np.abs(b._array))


def test_two_index_assign_two_index_transform():
    lf = DenseLinalgFactory(3)
    a = lf.create_two_index()
    e0 = lf.create_expansion()
    e1 = lf.create_expansion()
    a.randomize()
    a.symmetrize()
    e0.randomize()
    e1.randomize()
    b = a.new()
    b.assign_two_index_transform(a, e0)
    assert np.allclose(b._array, b._array.T)
    assert b.is_symmetric()
    assert np.allclose(b._array, np.dot(e0.coeffs.T, np.dot(a._array, e0.coeffs)))
    b.assign_two_index_transform(a, e0, e1)
    assert not b.is_symmetric()
    assert np.allclose(b._array, np.dot(e0.coeffs.T, np.dot(a._array, e1.coeffs)))


#
# DenseThreeIndex tests
#


def test_three_index_hdf5():
    lf = DenseLinalgFactory(5)
    a = lf.create_three_index()
    a.randomize()
    with h5.File('horton.matrix.test.test_dens.test_three_index_hdf5', driver='core', backing_store=False) as f:
        a.to_hdf5(f)
        b = DenseThreeIndex.from_hdf5(f)
        assert a == b


def test_three_index_copy_new_randomize_clear_assign():
    lf = DenseLinalgFactory(5)
    for args in (None,), (4,):
        a = lf.create_three_index(*args)
        b = a.copy()
        b.randomize()
        assert a != b
        c = b.copy()
        c.new.__check_init_args__(c, b)
        assert b == c
        d = c.new()
        assert a == d
        b.clear()
        assert a == b
        b.randomize()
        e = b.copy(begin0=1, end1=2, begin2=3)
        assert (e._array == b._array[1:,:2,3:]).all()


def test_three_index_permute_basis():
    lf = DenseLinalgFactory(5)
    for i in xrange(10):
        forth, back = get_forth_back(5)
        a = lf.create_three_index()
        a.randomize()
        b = a.copy()
        b.permute_basis(forth)
        assert a != b
        b.permute_basis(back)
        assert a == b


def test_three_index_change_basis_signs():
    lf = DenseLinalgFactory(5)
    for i in xrange(10):
        signs = get_signs(5)
        a = lf.create_three_index()
        a.randomize()
        b = a.copy()
        b.change_basis_signs(signs)
        assert a != b
        b.change_basis_signs(signs)
        assert a == b


def test_three_index_iadd():
    lf = DenseLinalgFactory(5)
    a = lf.create_three_index()
    a.randomize()
    b = lf.create_three_index()
    b.randomize()
    c = b.copy()
    factor = np.random.uniform(1, 2)
    c.iadd(a, factor)
    for i0 in xrange(5):
        for i1 in xrange(5):
            for i2 in xrange(5):
                assert factor*a.get_element(i0, i1, i2) + b.get_element(i0, i1, i2) == c.get_element(i0, i1, i2)


def test_three_index_iadd_slice():
    lf = DenseLinalgFactory(5)
    a = lf.create_three_index(9)
    a.randomize()
    b = lf.create_three_index()
    b.randomize()
    c = b.copy()
    factor = np.random.uniform(1, 2)
    c.iadd_slice(a, factor, 0, 5, 3, 8, 4, 9)
    for i0 in xrange(5):
        for i1 in xrange(5):
            for i2 in xrange(5):
                assert factor*a.get_element(i0, i1+3, i2+4) + b.get_element(i0, i1, i2) == c.get_element(i0, i1, i2)


def test_three_index_iscale():
    lf = DenseLinalgFactory()
    op = lf.create_three_index(3)
    op.randomize()
    tmp = op._array.copy()
    op.iscale(3.0)
    assert abs(op._array - 3*tmp).max() < 1e-10


def test_three_index_get_set():
    lf = DenseLinalgFactory()
    op = lf.create_three_index(3)
    op.set_element(0, 1, 2, 1.2)
    assert op.get_element(0, 1, 2) == 1.2


def test_three_index_contract_two_to_two():
    # Only test output for one case: abc,ab->ac
    lf = DenseLinalgFactory(3)
    three = lf.create_three_index()
    two = lf.create_two_index()
    three.randomize()
    two.randomize()
    out = three.contract_two_to_two('abc,ab->ac', two, factor=1.9)
    assert np.allclose(out._array, 1.9*np.einsum('abc,ab->ac', three._array, two._array))
    foo = three.contract_two_to_two('abc,ab->ac', two, out, factor=1.6)
    assert foo is out
    assert np.allclose(out._array, 1.6*np.einsum('abc,ab->ac', three._array, two._array))
    three.contract_two_to_two('abc,ab->ac', two, out, factor=1.6, clear=False)
    assert np.allclose(out._array, 3.2*np.einsum('abc,ab->ac', three._array, two._array))
    # Bindly test all other options
    three.contract_two_to_two('abc,ab->ca', two, factor=1.9)
    three.contract_two_to_two('abc,bc->ab', two, factor=1.9)
    three.contract_two_to_two('abc,bc->ba', two, factor=1.9)
    three.contract_two_to_two('abc,cb->ac', two, factor=1.9)


def test_three_index_iadd_expand_two_one():
    lf = DenseLinalgFactory(3)
    three = lf.create_three_index()
    two = lf.create_two_index()
    one = lf.create_one_index()
    three.randomize()
    two.randomize()
    one.randomize()
    orig = three.copy()
    three.iadd_expand_two_one('ab,c->cab', two, one, factor=0.7)
    assert np.allclose(three._array, orig._array + 0.7*np.einsum('ab,c->cab', two._array, one._array))
    orig = three.copy()
    three.iadd_expand_two_one('ab,c->acb', two, one, factor=0.7)
    assert np.allclose(three._array, orig._array + 0.7*np.einsum('ab,c->acb', two._array, one._array))


def test_three_index_iadd_expand_two_two():
    # Only test output for two cases: 'ac,bc->abc', 'cb,ac->acb'
    lf = DenseLinalgFactory(3)
    three = lf.create_three_index()
    two1 = lf.create_two_index()
    two2 = lf.create_two_index()
    three.randomize()
    two1.randomize()
    two2.randomize()
    orig = three.copy()
    three.iadd_expand_two_two('ac,bc->abc', two1, two2, factor=0.7)
    assert np.allclose(three._array, orig._array + 0.7*np.einsum('ac,bc->abc', two1._array, two2._array))
    # blind test for remaining contractions
    others = ['ac,bc->abc', 'ab,bc->abc', 'ab,ac->acb', 'cb,ac->acb',
              'ac,ab->abc', 'ab,ac->abc']
    for select in others:
        three.iadd_expand_two_two(select, two1, two2, factor=0.7)
    # with ranges
    two1 = lf.create_two_index(5)
    two2 = lf.create_two_index(7)
    two1.randomize()
    two2.randomize()
    orig = three.copy()
    three.iadd_expand_two_two('cb,ac->acb', two1, two2, factor=0.7, end0=3, begin1=3, end1=6, begin2=2, end2=5, begin3=2)
    assert np.allclose(three._array, orig._array + 0.7*np.einsum('cb,ac->acb', two1._array[2:5,2:], two2._array[:3,3:6]))
    # blind test for remaining contractions
    others = ['ac,bc->abc', 'ab,bc->abc', 'ab,ac->acb', 'cb,ac->acb',
              'ac,ab->abc', 'ab,ac->abc']
    for select in others:
        three.iadd_expand_two_two(select, two1, two2, factor=0.7, end0=3, begin1=3, end1=6, begin2=2, end2=5, begin3=2)


def test_three_index_contract_two_to_three():
    # Test only for one case of subscripts
    lf = DenseLinalgFactory(3)
    a = lf.create_three_index()
    b = lf.create_three_index()
    two = lf.create_two_index()
    a.randomize()
    b.randomize()
    two.randomize()
    orig = b.copy()
    a.contract_two_to_three('abc,ad->cdb', two, b, factor=0.7, clear=False)
    assert np.allclose(b._array, orig._array + 0.7*np.einsum('abc,ad->cdb', a._array, two._array))
    out = a.contract_two_to_three('abc,ad->cdb', two, b, factor=0.7, clear=True)
    assert out is b
    assert np.allclose(b._array, 0.7*np.einsum('abc,ad->cdb', a._array, two._array))
    # try other cases blindly
    others = ['abc,ad->cbd', 'abc,da->dbc', 'abc,da->bdc']
    for select in others:
        a.contract_two_to_three(select, two, b, factor=0.7)
    # with ranges
    a = lf.create_three_index(9)
    a.randomize()
    out = a.contract_two_to_three('abc,ad->cdb', two, b, factor=0.7, clear=True, end0=3, begin1=2, end1=5, begin2=6, end2=9)
    assert out is b
    assert np.allclose(b._array, 0.7*np.einsum('abc,ad->cdb', a._array[:3, 2:5, 6:9], two._array))
    out = a.contract_two_to_three('abc,ad->cdb', two, b, factor=0.7, clear=False, end0=3, begin1=2, end1=5, begin2=6, end2=9)
    assert np.allclose(b._array, 1.4*np.einsum('abc,ad->cdb', a._array[:3, 2:5, 6:9], two._array))


def test_three_index_contract_to_two():
    # Test only for one case of subscripts
    lf = DenseLinalgFactory(3)
    a = lf.create_three_index()
    b = lf.create_two_index()
    a.randomize()
    b.randomize()
    out = a.contract_to_two('abc->ac', b, factor=0.7, clear=True)
    assert out is b
    assert np.allclose(b._array, 0.7*np.einsum('abc->ac', a._array))
    orig = b.copy()
    out = a.contract_to_two('abc->ac', b, factor=0.7, clear=False)
    assert out is b
    assert np.allclose(b._array, orig._array + 0.7*np.einsum('abc->ac', a._array))
    # with ranges
    a = lf.create_three_index(9)
    a.randomize()
    out = a.contract_to_two('abc->ac', b, factor=0.7, clear=True, end0=3, begin1=2, end1=5, begin2=6, end2=9)
    assert np.allclose(b._array, 0.7*np.einsum('abc->ac', a._array[:3, 2:5, 6:9]))
    orig = b.copy()
    out = a.contract_to_two('abc->ac', b, factor=0.7, clear=False, end0=3, begin1=2, end1=5, begin2=6, end2=9)
    assert out is b
    assert np.allclose(b._array, orig._array + 0.7*np.einsum('abc->ac', a._array[:3, 2:5, 6:9]))


#
# DenseFourIndex tests
#

# FIXME: extend tests for different symmetries

def test_four_index_hdf5():
    lf = DenseLinalgFactory(5)
    a = lf.create_four_index()
    a.randomize()
    with h5.File('horton.matrix.test.test_dens.test_four_index_hdf5', driver='core', backing_store=False) as f:
        a.to_hdf5(f)
        b = DenseFourIndex.from_hdf5(f)
        assert a == b


def test_four_index_copy_new_randomize_clear_assign():
    lf = DenseLinalgFactory(5)
    for args in (None,), (4,):
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
        b.randomize()
        e = b.copy(0, 1, 2, 4, 2)
        assert (e._array == b._array[:1,2:4,2:]).all()


def test_four_index_permute_basis():
    lf = DenseLinalgFactory(5)
    for i in xrange(10):
        forth, back = get_forth_back(5)
        a = lf.create_four_index()
        a.randomize()
        b = a.copy()
        b.permute_basis(forth)
        assert a != b
        b.permute_basis(back)
        assert a == b


def test_four_index_change_basis_signs():
    lf = DenseLinalgFactory(5)
    for i in xrange(10):
        signs = get_signs(5)
        a = lf.create_four_index()
        a.randomize()
        b = a.copy()
        b.change_basis_signs(signs)
        assert a != b
        b.change_basis_signs(signs)
        assert a == b


def test_four_index_iadd():
    lf = DenseLinalgFactory(5)
    a = lf.create_four_index()
    a.randomize()
    b = lf.create_four_index()
    b.randomize()
    c = b.copy()
    factor = np.random.uniform(1, 2)
    c.iadd(a, factor)
    for i0 in xrange(5):
        for i1 in xrange(5):
            for i2 in xrange(5):
                for i3 in xrange(5):
                    assert factor*a.get_element(i0, i1, i2, i3) + b.get_element(i0, i1, i2, i3) == c.get_element(i0, i1, i2, i3)


def test_four_index_iscale():
    lf = DenseLinalgFactory()
    op = lf.create_four_index(3)
    op.randomize()
    tmp = op._array.copy()
    op.iscale(3.0)
    assert abs(op._array - 3*tmp).max() < 1e-10


def test_four_index_get_set():
    lf = DenseLinalgFactory()
    op = lf.create_four_index(4)
    op.set_element(0, 1, 2, 3, 1.2)
    assert op.get_element(0, 1, 2, 3) == 1.2
    assert op.get_element(2, 1, 0, 3) == 1.2
    assert op.get_element(0, 3, 2, 1) == 1.2
    assert op.get_element(2, 3, 0, 1) == 1.2
    assert op.get_element(1, 0, 3, 2) == 1.2
    assert op.get_element(3, 0, 1, 2) == 1.2
    assert op.get_element(1, 2, 3, 0) == 1.2
    assert op.get_element(3, 2, 1, 0) == 1.2


def test_four_index_is_symmetric():
    lf = DenseLinalgFactory(4)
    op = lf.create_four_index()
    op.set_element(0, 1, 2, 3, 1.234)
    assert op.is_symmetric(8)
    assert op.is_symmetric(4)
    assert op.is_symmetric(2)
    assert op.is_symmetric(1)
    op.set_element(0, 1, 2, 3, 1.0, symmetry=4)
    assert not op.is_symmetric(8)
    assert not op.is_symmetric(2)
    assert op.is_symmetric(4)
    assert op.is_symmetric(1)
    op.set_element(0, 1, 2, 3, 1.234)
    op.set_element(0, 1, 2, 3, 0.5, symmetry=2)
    assert not op.is_symmetric(8)
    assert not op.is_symmetric(4)
    assert op.is_symmetric(2)
    assert op.is_symmetric(1)
    op.set_element(0, 1, 2, 3, 0.3, symmetry=1)
    assert not op.is_symmetric(8)
    assert not op.is_symmetric(4)
    assert not op.is_symmetric(2)
    assert op.is_symmetric(1)


def test_four_index_symmetrize():
    lf = DenseLinalgFactory(8)
    op = lf.create_four_index()
    for symmetry in 1, 2, 4, 8:
        op.randomize()
        op.symmetrize(symmetry)
        assert op.is_symmetric(symmetry, 0, 0)


def test_four_index_symmetrize_order_of_operations():
    lf = DenseLinalgFactory(8)
    op = lf.create_four_index()
    for symmetry in 1, 2, 4, 8:
        # ugly hack to have matrix elements with very different order of
        # magnitudes
        op._array[:] = 10**np.random.uniform(-20,20, (8,8,8,8))
        op.symmetrize(symmetry)
        assert op.is_symmetric(symmetry, 0, 0)


def test_four_index_itranspose():
    lf = DenseLinalgFactory(8)
    for i in xrange(10):
        op = lf.create_four_index()
        op.randomize()
        i0, i1, i2, i3 = np.random.randint(0, 4, 4)
        x = op.get_element(i0, i1, i2, i3)
        op.itranspose()
        assert op.get_element(i1, i0, i3, i2) == x


def test_four_index_sum():
    # Blind test
    lf = DenseLinalgFactory(4)
    op = lf.create_four_index()
    op.randomize()
    op.sum()


def test_four_index_iadd_exchange():
    # Blind test
    lf = DenseLinalgFactory(4)
    op = lf.create_four_index()
    op.randomize()
    op.symmetrize()
    op.iadd_exchange()
    op.is_symmetric()


def test_four_index_slice_to_two():
    # test in detail for aabb->ab
    lf = DenseLinalgFactory(6)
    four = lf.create_four_index()
    four.randomize()
    two = four.slice_to_two('aabb->ab', factor=1.3)
    assert np.allclose(two._array, 1.3*np.einsum('aabb->ab', four._array))
    foo = four.slice_to_two('aabb->ab', two, factor=1.4)
    assert foo is two
    assert np.allclose(two._array, 1.4*np.einsum('aabb->ab', four._array))
    four.slice_to_two('aabb->ab', two, factor=1.4, clear=False)
    assert np.allclose(two._array, 2.8*np.einsum('aabb->ab', four._array))
    # Blind test on all other cases
    four.slice_to_two('abab->ab', factor=1.3)
    four.slice_to_two('abba->ab', factor=1.3)
    # with ranges
    two = four.slice_to_two('aabb->ab', factor=1.3, end0=3, end1=3, begin2=2, begin3=2)
    assert np.allclose(two._array, 1.3*np.einsum('aabb->ab', four._array[:3,:3,2:,2:]))
    foo = four.slice_to_two('aabb->ab', two, factor=1.4, end0=3, end1=3, begin2=2, begin3=2)
    assert foo is two
    assert np.allclose(two._array, 1.4*np.einsum('aabb->ab', four._array[:3,:3,2:,2:]))
    four.slice_to_two('aabb->ab', two, factor=1.4, clear=False, end0=3, end1=3, begin2=2, begin3=2)
    assert np.allclose(two._array, 2.8*np.einsum('aabb->ab', four._array[:3,:3,2:,2:]))
    # Blind test on all other cases
    four.slice_to_two('abab->ab', factor=1.3, end0=3, begin1=2, end2=3, begin3=2)
    four.slice_to_two('abba->ab', factor=1.3, end0=3, begin1=2, begin2=2, end3=3)


def test_four_index_slice_to_three():
    # test in detail for aabb->ab
    lf = DenseLinalgFactory(4)
    four = lf.create_four_index()
    four.randomize()
    three = four.slice_to_three('abcc->bac', factor=1.3)
    assert np.allclose(three._array, 1.3*np.einsum('abcc->bac', four._array))
    foo = four.slice_to_three('abcc->bac', three, factor=1.4)
    assert foo is three
    assert np.allclose(three._array, 1.4*np.einsum('abcc->bac', four._array))
    four.slice_to_three('abcc->bac', three, factor=1.4, clear=False)
    assert np.allclose(three._array, 2.8*np.einsum('abcc->bac', four._array))
    # Blind test on all other cases
    four.slice_to_three('abcc->abc', factor=1.3)
    four.slice_to_three('abcb->abc', factor=1.3)
    four.slice_to_three('abbc->abc', factor=1.3)


def test_four_index_slice_to_four():
    # test in detail for abcd->abcd'
    lf = DenseLinalgFactory(6)
    four = lf.create_four_index()
    four.randomize()
    four2 = four.slice_to_four('abcd->abcd', factor=1.3)
    assert np.allclose(four2._array, 1.3*np.einsum('abcd->abcd', four._array))
    foo = four.slice_to_four('abcd->abcd', four2, factor=1.4)
    assert foo is four2
    assert np.allclose(four2._array, 1.4*np.einsum('abcd->abcd', four._array))
    four.slice_to_four('abcd->abcd', four2, factor=1.4, clear=False)
    assert np.allclose(four2._array, 2.8*np.einsum('abcd->abcd', four._array))
    # Blind test on all other cases
    four.slice_to_four('abcd->acbd', factor=1.3)
    four.slice_to_four('abcd->cadb', factor=1.3)
    # with ranges
    four2 = four.slice_to_four('abcd->abcd', factor=1.3, end0=3, begin1=2, begin2=2, end2=5)
    assert np.allclose(four2._array, 1.3*np.einsum('abcd->abcd', four._array[:3,2:,2:5,:]))
    foo = four.slice_to_four('abcd->abcd', four2, factor=1.4, end0=3, begin1=2, begin2=2, end2=5)
    assert foo is four2
    assert np.allclose(four2._array, 1.4*np.einsum('abcd->abcd', four._array[:3,2:,2:5,:]))
    four.slice_to_four('abcd->abcd', four2, factor=1.4, clear=False, end0=3, begin1=2, begin2=2, end2=5)
    assert np.allclose(four2._array, 2.8*np.einsum('abcd->abcd', four._array[:3,2:,2:5,:]))
    # Blind test on all other cases
    four.slice_to_four('abcd->acbd', factor=1.3, end0=3, begin1=2, begin2=2, end2=5)
    four.slice_to_four('abcd->cadb', factor=1.3, end0=3, begin1=2, begin2=2, end2=5)


old_numpy_warning = """
   Your NumPy version is too old. Please install version 1.9.1 or newer. Consult
   to the "Download and install" documentation for more details.
"""


def test_four_index_contract_two():
    # test in detail for aabb,ab
    lf = DenseLinalgFactory(6)
    inp = lf.create_four_index()
    two = lf.create_two_index()
    inp.randomize()
    two.randomize()
    try:
        out = inp.contract_two('aabb,ab', two, factor=1.3)
    except ValueError:
        if np.__version__ < '1.9.1':
            print old_numpy_warning
        raise
    assert out == 1.3*np.einsum('aabb,ab', inp._array, two._array)
    # with ranges
    two = lf.create_two_index(3)
    two.randomize()
    out = inp.contract_two('aabb,ab', two, factor=1.3, end0=3, end1=3, begin2=2, begin3=2, end2=5, end3=5)
    assert out == 1.3*np.einsum('aabb,ab', inp._array[:3,:3,2:5,2:5], two._array)


def test_four_index_contract_four():
    # test in detail for abcd,abcd
    lf = DenseLinalgFactory(3)
    inp = lf.create_four_index()
    four = lf.create_four_index()
    inp.randomize()
    four.randomize()
    out = inp.contract_four('abcd,abcd', four, factor=1.3)
    assert out == 1.3*np.einsum('abcd,abcd', inp._array, four._array)
    # Blind test on all other cases
    others = ['abcd,abcd', 'abcd,adcb', 'abab,abab', 'abad,abad', 'abdb,abdb',
              'abcd,acbd', 'abcd,acdb']
    for subscripts in others:
        inp.contract_four(subscripts, four, factor=1.7)
    # with ranges
    four = lf.create_four_index(5)
    four.randomize()
    out = inp.contract_four('abcd,abcd', four, factor=1.3, end0=3, begin1=2, begin2=1, end2=4, begin3=2)
    assert out == 1.3*np.einsum('abcd,abcd', inp._array, four._array[:3,2:,1:4,2:])
    # Blind test on all other cases
    others = ['abcd,abcd', 'abcd,adcb', 'abab,abab', 'abad,abad', 'abdb,abdb',
              'abcd,acbd', 'abcd,acdb']
    for subscripts in others:
        inp.contract_four(subscripts, four, factor=1.7, end0=3, begin1=2, begin2=1, end2=4, begin3=2)


def test_four_index_contract_to_two():
    # test in detail for 'abcb->ac', 'abbc->ac'
    lf = DenseLinalgFactory(5)
    inp = lf.create_four_index()
    two = lf.create_two_index()
    inp.randomize()
    two.randomize()
    out = inp.contract_to_two('abcb->ac', two, factor=1.3)
    assert np.allclose(out._array, 1.3*np.einsum('abcb->ac', inp._array))
    foo = inp.contract_to_two('abcb->ac', out, factor=1.4)
    assert foo is out
    assert np.allclose(out._array, 1.4*np.einsum('abcb->ac', inp._array))
    inp.contract_to_two('abcb->ac', out, factor=1.4, clear=False)
    assert np.allclose(out._array, 2.8*np.einsum('abcb->ac', inp._array))
    # Blind test on all other cases
    others = ['abbc->ac']
    for subscripts in others:
        inp.contract_to_two(subscripts, two, factor=1.7)
    # with ranges
    two = lf.create_two_index(3)
    inp.randomize()
    out = inp.contract_to_two('abcb->ac', two, factor=1.3, end0=3, begin1=2, begin2=1, end2=4, begin3=2)
    assert np.allclose(out._array, 1.3*np.einsum('abcb->ac', inp._array[:3,2:,1:4,2:]))
    foo = inp.contract_to_two('abcb->ac', out, factor=1.4, end0=3, begin1=2, begin2=1, end2=4, begin3=2)
    assert foo is out
    assert np.allclose(out._array, 1.4*np.einsum('abcb->ac', inp._array[:3,2:,1:4,2:]))
    inp.contract_to_two('abcb->ac', out, factor=1.4, clear=False, end0=3, begin1=2, begin2=1, end2=4, begin3=2)
    assert np.allclose(out._array, 2.8*np.einsum('abcb->ac', inp._array[:3,2:,1:4,2:]))
    # Blind test on all other cases
    others = ['abbc->ac']
    for subscripts in others:
        inp.contract_to_two(subscripts, two, factor=1.7, end0=3, begin1=2, begin2=1, end2=4, begin3=2)


def test_four_index_contract_two_to_four():
    # test in detail for abcd,cd->acbd
    lf = DenseLinalgFactory(6)
    inp = lf.create_four_index()
    two = lf.create_two_index()
    inp.randomize()
    two.randomize()
    out = inp.contract_two_to_four('abcd,cd->acbd', two, factor=1.3)
    assert np.allclose(out._array, 1.3*np.einsum('abcd,cd->acbd', inp._array, two._array))
    foo = inp.contract_two_to_four('abcd,cd->acbd', two, out, factor=1.4)
    assert foo is out
    assert np.allclose(out._array, 1.4*np.einsum('abcd,cd->acbd', inp._array, two._array))
    inp.contract_two_to_four('abcd,cd->acbd', two, out, factor=1.4, clear=False)
    assert np.allclose(out._array, 2.8*np.einsum('abcd,cd->acbd', inp._array, two._array))
    # Blind test on all other cases
    others = ['abcd,cd->acdb', 'abcd,cb->acdb', 'abcd,cb->acbd', 'abcd,ab->acbd',
            'abcd,ab->acdb', 'abcd,ad->acbd', 'abcd,ad->acdb', 'abcd,ad->abcd',
            'abcd,ad->abdc', 'abcd,bd->abcd', 'abcd,bd->abdc', 'abcd,bc->abdc',
            'abcd,bc->abcd', 'abcd,ac->abcd', 'abcd,ac->abdc', 'abcd,bd->cabd',
            'abcd,bc->dabc', 'abcd,ca->cabd', 'abcd,da->dabc', 'abcd,dc->dabc',
            'abcd,ba->dabc', 'abcd,ac->acbd', 'abcd,cd->abcd', 'abcc,ad->abcd',
            'aabc,dc->adbc', 'aabc,db->adbc', 'abcc,ad->abcd', 'abcc,bd->abcd',
            'abcd,bc->acbd', 'abcd,eb->aecd', 'abcd,eb->cdae', 'abcd,ed->ceab',
            'abcd,ed->abce', 'abcd,ae->cdeb', 'abcd,ae->ebcd', 'abcd,ce->edab',
            'abcd,ce->abed', 'abcd,ab->abcd', 'abcd,cb->abcd', 'abcd,ec->eadb',
            'abcd,ec->dbea', 'abcd,ae->cedb', 'abcd,ae->dbce', 'abcd,ec->ebad',
            'abcd,ed->ebac', 'abcd,ec->adeb', 'abcd,ed->aceb', 'abcd,ae->debc',
            'abcd,ae->cebd', 'abcd,ae->bcde', 'abcd,ae->bdce']
    for subscripts in others:
        inp.contract_two_to_four(subscripts, two, factor=1.7)
    # with ranges
    out = inp.contract_two_to_four('abcd,cd->acbd', two, factor=1.3,
        begin0=0, end0=3, begin1=1, end1=4, begin2=2, end2=5,
        begin3=3, end3=6, begin4=2, end4=5, begin5=3, end5=6)
    assert np.allclose(out._array, 1.3*np.einsum('abcd,cd->acbd', inp._array[:3,1:4,2:5,3:6], two._array[2:5,3:6]))
    foo = inp.contract_two_to_four('abcd,cd->acbd', two, out, factor=1.4,
        begin0=0, end0=3, begin1=1, end1=4, begin2=2, end2=5,
        begin3=3, end3=6, begin4=2, end4=5, begin5=3, end5=6)
    assert foo is out
    assert np.allclose(out._array, 1.4*np.einsum('abcd,cd->acbd', inp._array[:3,1:4,2:5,3:6], two._array[2:5,3:6]))
    inp.contract_two_to_four('abcd,cd->acbd', two, out, factor=1.4, clear=False,
        begin0=0, end0=3, begin1=1, end1=4, begin2=2, end2=5,
        begin3=3, end3=6, begin4=2, end4=5, begin5=3, end5=6)
    assert np.allclose(out._array, 2.8*np.einsum('abcd,cd->acbd', inp._array[:3,1:4,2:5,3:6], two._array[2:5,3:6]))
    for subscripts in others:
        inp.contract_two_to_four(subscripts, two, factor=1.7,
            begin0=0, end0=3, begin1=1, end1=4, begin2=2, end2=5,
            begin3=3, end3=6, begin4=2, end4=5, begin5=3, end5=6)


def test_four_index_contract_two_to_two():
    # test in detail for abcd,bd->ac
    lf = DenseLinalgFactory(6)
    inp = lf.create_four_index()
    two = lf.create_two_index()
    inp.randomize()
    inp.symmetrize()
    two.randomize()
    two.symmetrize()
    out = inp.contract_two_to_two('abcd,bd->ac', two, factor=1.3)
    assert np.allclose(out._array, 1.3*np.einsum('abcd,bd->ac', inp._array, two._array))
    foo = inp.contract_two_to_two('abcd,bd->ac', two, out, factor=1.4)
    assert foo is out
    assert np.allclose(out._array, 1.4*np.einsum('abcd,bd->ac', inp._array, two._array))
    inp.contract_two_to_two('abcd,bd->ac', two, out, factor=1.4, clear=False)
    assert np.allclose(out._array, 2.8*np.einsum('abcd,bd->ac', inp._array, two._array))
    # Blind test on all other cases
    others = ['abcd,cb->ad',
              'aabb,cb->ac', 'aabb,cb->ca', 'abcc,bc->ab', 'aabc,ab->bc',
              'aabc,ac->bc', 'abcc,ac->ab', 'abcb,cb->ac', 'abcb,ab->ac',
              'abcc,bc->ba', 'abcc,bc->ab', 'abcd,ac->db', 'abcd,ad->cb',
              'abcd,ac->bd', 'abcd,ad->bc', 'abcd,ab->cd']
    for subscripts in others:
        try:
            inp.contract_two_to_two(subscripts, two, factor=1.7)
        except ValueError:
            if np.__version__ < '1.9.1':
                print old_numpy_warning
            raise
    # with ranges
    out = inp.contract_two_to_two('abcd,cb->ad', two, factor=1.3,
        begin0=0, end0=3, begin1=1, end1=4, begin2=2, end2=5,
        begin3=3, end3=6, begin4=2, end4=5, begin5=3, end5=6)
    assert np.allclose(out._array, 1.3*np.einsum('abcd,cb->ad', inp._array[:3,1:4,2:5,3:6], two._array[2:5,3:6]))
    foo = inp.contract_two_to_two('abcd,cb->ad', two, out, factor=1.4,
        begin0=0, end0=3, begin1=1, end1=4, begin2=2, end2=5,
        begin3=3, end3=6, begin4=2, end4=5, begin5=3, end5=6)
    assert foo is out
    assert np.allclose(out._array, 1.4*np.einsum('abcd,cb->ad', inp._array[:3,1:4,2:5,3:6], two._array[2:5,3:6]))
    inp.contract_two_to_two('abcd,cb->ad', two, out, factor=1.4, clear=False,
        begin0=0, end0=3, begin1=1, end1=4, begin2=2, end2=5,
        begin3=3, end3=6, begin4=2, end4=5, begin5=3, end5=6)
    assert np.allclose(out._array, 2.8*np.einsum('abcd,cb->ad', inp._array[:3,1:4,2:5,3:6], two._array[2:5,3:6]))
    # Blind test on all other cases
    for subscripts in others:
        inp.contract_two_to_two('abcd,cb->ad', two, out, factor=1.7,
            begin0=0, end0=3, begin1=1, end1=4, begin2=2, end2=5,
            begin3=3, end3=6, begin4=2, end4=5, begin5=3, end5=6)


def test_four_index_contract_two_to_three():
    # test in detail for 'aabc,ad->bcd'
    lf = DenseLinalgFactory(6)
    inp = lf.create_four_index()
    two = lf.create_two_index()
    inp.randomize()
    two.randomize()
    out = inp.contract_two_to_three('aabc,ad->bcd', two, factor=1.3)
    assert np.allclose(out._array, 1.3*np.einsum('aabc,ad->bcd', inp._array, two._array))
    foo = inp.contract_two_to_three('aabc,ad->bcd', two, out, factor=1.4)
    assert foo is out
    assert np.allclose(out._array, 1.4*np.einsum('aabc,ad->bcd', inp._array, two._array))
    inp.contract_two_to_three('aabc,ad->bcd', two, out, factor=1.4, clear=False)
    assert np.allclose(out._array, 2.8*np.einsum('aabc,ad->bcd', inp._array, two._array))
    # Blind test on all other cases
    others = [ 'abcd,ac->bdc',
               'abcd,ad->bcd', 'abcd,bc->abd', 'abcd,ac->abd', 'abcc,dc->dab',
               'abcd,ac->bcd', 'abcc,cd->dab', 'abcc,dc->abd', 'abcb,db->adc',
               'abcc,dc->adb', 'abcb,db->dac',
               'abbc,db->dac', 'abcc,cd->abd']
    for subscripts in others:
        inp.contract_two_to_three(subscripts, two, factor=1.7)
    # with ranges
    two = lf.create_two_index(3)
    two.randomize()
    out = inp.contract_two_to_three('aabc,ad->bcd', two, factor=1.3,
            begin0=0, end0=3, begin1=1, end1=4, begin2=2, end2=5,
            begin3=3, end3=6)
    assert np.allclose(out._array, 1.3*np.einsum('aabc,ad->bcd', inp._array[:3,1:4,2:5,3:6], two._array))
    foo = inp.contract_two_to_three('aabc,ad->bcd', two, out, factor=1.4,
            begin0=0, end0=3, begin1=1, end1=4, begin2=2, end2=5,
            begin3=3, end3=6)
    assert foo is out
    assert np.allclose(out._array, 1.4*np.einsum('aabc,ad->bcd', inp._array[:3,1:4,2:5,3:6], two._array))
    inp.contract_two_to_three('aabc,ad->bcd', two, out, factor=1.4, clear=False,
            begin0=0, end0=3, begin1=1, end1=4, begin2=2, end2=5,
            begin3=3, end3=6)
    assert np.allclose(out._array, 2.8*np.einsum('aabc,ad->bcd', inp._array[:3,1:4,2:5,3:6], two._array))
    # Blind test on all other cases
    for subscripts in others:
        inp.contract_two_to_three(subscripts, two, factor=1.7,
            begin0=0, end0=3, begin1=1, end1=4, begin2=2, end2=5,
            begin3=3, end3=6)


def test_four_index_contract_three_to_three():
    # test in detail for 'abcd,ace->ebd',
    lf = DenseLinalgFactory(3)
    inp = lf.create_four_index()
    three = lf.create_three_index()
    inp.randomize()
    three.randomize()
    out = inp.contract_three_to_three('abcd,ace->ebd', three, factor=1.3)
    assert np.allclose(out._array, 1.3*np.einsum('abcd,ace->ebd', inp._array, three._array))
    foo = inp.contract_three_to_three('abcd,ace->ebd', three, out, factor=1.4)
    assert foo is out
    assert np.allclose(out._array, 1.4*np.einsum('abcd,ace->ebd', inp._array, three._array))
    inp.contract_three_to_three('abcd,ace->ebd', three, out, factor=1.4, clear=False)
    assert np.allclose(out._array, 2.8*np.einsum('abcd,ace->ebd', inp._array, three._array))
    # Blind test on all other cases
    others = [ 'abcd,ebd->ace']
    for subscripts in others:
        inp.contract_three_to_three(subscripts, three, factor=1.7)
    # with ranges
    three = lf.create_three_index(6)
    three.randomize()
    out = inp.contract_three_to_three('abcd,ace->ebd', three, factor=1.3,
            begin0=0, end0=3, begin1=1, end1=4, begin2=2, end2=5)
    assert np.allclose(out._array, 1.3*np.einsum('abcd,ace->ebd', inp._array, three._array[:3,1:4,2:5]))
    foo = inp.contract_three_to_three('abcd,ace->ebd', three, out, factor=1.4,
            begin0=0, end0=3, begin1=1, end1=4, begin2=2, end2=5)
    assert foo is out
    assert np.allclose(out._array, 1.4*np.einsum('abcd,ace->ebd', inp._array, three._array[:3,1:4,2:5]))
    inp.contract_three_to_three('abcd,ace->ebd', three, out, factor=1.4, clear=False,
            begin0=0, end0=3, begin1=1, end1=4, begin2=2, end2=5)
    assert np.allclose(out._array, 2.8*np.einsum('abcd,ace->ebd', inp._array, three._array[:3,1:4,2:5]))
    for subscripts in others:
        inp.contract_three_to_three(subscripts, three, factor=1.7,
            begin0=0, end0=3, begin1=1, end1=4, begin2=2, end2=5)


def test_four_index_contract_four_to_two():
    # test in detail for 'abcd,aced->be',
    lf = DenseLinalgFactory(3)
    inp = lf.create_four_index()
    four = lf.create_four_index()
    inp.randomize()
    four.randomize()
    out = inp.contract_four_to_two('abcd,aced->be', four, factor=1.3)
    assert np.allclose(out._array, 1.3*np.einsum('abcd,aced->be', inp._array, four._array))
    foo = inp.contract_four_to_two('abcd,aced->be', four, out, factor=1.4)
    assert foo is out
    assert np.allclose(out._array, 1.4*np.einsum('abcd,aced->be', inp._array, four._array))
    inp.contract_four_to_two('abcd,aced->be', four, out, factor=1.4, clear=False)
    assert np.allclose(out._array, 2.8*np.einsum('abcd,aced->be', inp._array, four._array))
    # Blind test on all other cases
    others = ['abcd,acde->be', 'abcd,aced->eb', 'abcd,acde->eb',
              'abcd,ecbd->ae', 'abcd,ecdb->ae', 'abcd,ecdb->ea',
              'abcd,ecbd->ea', 'abcd,acbe->ed', 'abcd,aceb->ed',
              'abcd,aebd->ce', 'abcd,aedb->ce']
    for subscripts in others:
        inp.contract_four_to_two(subscripts, four, factor=1.7)
    # with ranges
    four = lf.create_four_index(6)
    four.randomize()
    out = inp.contract_four_to_two('abcd,aced->be', four, factor=1.3,
            begin0=0, end0=3, begin1=1, end1=4, begin2=2, end2=5,
            begin3=3, end3=6)
    assert np.allclose(out._array, 1.3*np.einsum('abcd,aced->be', inp._array, four._array[:3,1:4,2:5,3:6]))
    foo = inp.contract_four_to_two('abcd,aced->be', four, out, factor=1.4,
            begin0=0, end0=3, begin1=1, end1=4, begin2=2, end2=5,
            begin3=3, end3=6)
    assert foo is out
    assert np.allclose(out._array, 1.4*np.einsum('abcd,aced->be', inp._array, four._array[:3,1:4,2:5,3:6]))
    inp.contract_four_to_two('abcd,aced->be', four, out, factor=1.4, clear=False,
            begin0=0, end0=3, begin1=1, end1=4, begin2=2, end2=5,
            begin3=3, end3=6)
    assert np.allclose(out._array, 2.8*np.einsum('abcd,aced->be', inp._array, four._array[:3,1:4,2:5,3:6]))
    # Blind test on all other cases
    for subscripts in others:
        inp.contract_four_to_two(subscripts, four, factor=1.7,
            begin0=0, end0=3, begin1=1, end1=4, begin2=2, end2=5,
            begin3=3, end3=6)


def test_four_index_contract_four_to_four():
    # test in detail for 'abcd,cedf->abfe',
    lf = DenseLinalgFactory(3)
    inp = lf.create_four_index()
    four = lf.create_four_index()
    inp.randomize()
    four.randomize()
    out = inp.contract_four_to_four('abcd,cedf->abfe', four, factor=1.3)
    assert np.allclose(out._array, 1.3*np.einsum('abcd,cedf->abfe', inp._array, four._array))
    foo = inp.contract_four_to_four('abcd,cedf->abfe', four, out, factor=1.4)
    assert foo is out
    assert np.allclose(out._array, 1.4*np.einsum('abcd,cedf->abfe', inp._array, four._array))
    inp.contract_four_to_four('abcd,cedf->abfe', four, out, factor=1.4, clear=False)
    assert np.allclose(out._array, 2.8*np.einsum('abcd,cedf->abfe', inp._array, four._array))
    # Blind test on all other cases
    others = ['abcd,cefd->abfe', 'abcd,cedf->feab', 'abcd,cefd->feab',
              'abcd,eadf->fbce', 'abcd,eadf->cefb', 'abcd,aedf->cbfe',
              'abcd,aedf->fecb', 'abcd,acef->ebfd', 'abcd,bdef->aecf',
              'abcd,cedf->abef', 'abcd,cefd->abef', 'abcd,cedf->efab',
              'abcd,cefd->efab', 'abcd,eadf->ebcf', 'abcd,eadf->cfeb',
              'abcd,eadf->cbef', 'abcd,eadf->efcb', 'abcd,eafd->cbef',
              'abcd,eafd->efcb']
    for subscripts in others:
        inp.contract_four_to_four(subscripts, four, factor=1.7)
    # with ranges
    four = lf.create_four_index(6)
    four.randomize()
    out = inp.contract_four_to_four('abcd,cedf->abfe', four, factor=1.3,
            begin0=0, end0=3, begin1=1, end1=4, begin2=2, end2=5,
            begin3=3, end3=6)
    assert np.allclose(out._array, 1.3*np.einsum('abcd,cedf->abfe', inp._array, four._array[:3,1:4,2:5,3:6]))
    foo = inp.contract_four_to_four('abcd,cedf->abfe', four, out, factor=1.4,
            begin0=0, end0=3, begin1=1, end1=4, begin2=2, end2=5,
            begin3=3, end3=6)
    assert foo is out
    assert np.allclose(out._array, 1.4*np.einsum('abcd,cedf->abfe', inp._array, four._array[:3,1:4,2:5,3:6]))
    inp.contract_four_to_four('abcd,cedf->abfe', four, out, factor=1.4, clear=False,
            begin0=0, end0=3, begin1=1, end1=4, begin2=2, end2=5,
            begin3=3, end3=6)
    assert np.allclose(out._array, 2.8*np.einsum('abcd,cedf->abfe', inp._array, four._array[:3,1:4,2:5,3:6]))
    # Blind test on all other cases
    for subscripts in others:
        inp.contract_four_to_four(subscripts, four, factor=1.7,
            begin0=0, end0=3, begin1=1, end1=4, begin2=2, end2=5,
            begin3=3, end3=6)


def test_four_index_iadd_expand_two_to_four():
    # test in detail for 0-2
    lf = DenseLinalgFactory(4)
    inp = lf.create_four_index()
    two = lf.create_two_index()
    inp.randomize()
    out = inp.copy()
    two.randomize()
    inp.iadd_expand_two_to_four('0-2', two, factor=1.3)
    for i in range(4):
        for j in range(4):
            for k in range(4):
                out._array[j,i,k,i] += two._array[j,k]*1.3
    assert np.allclose(out._array, inp._array)
    # Blind test on all other cases
    others = ['1-3', '0-3', '1-2', 'diag']
    for subscripts in others:
        inp.iadd_expand_two_to_four(subscripts, two, factor=1.7)
    # with ranges
    two = lf.create_two_index(8)
    two.randomize()
    out = inp.copy()
    inp.iadd_expand_two_to_four('0-2', two, factor=1.3, begin0=3, end0=7, end1=4)
    for i in range(4):
        for j in range(4):
            for k in range(4):
                out._array[j,i,k,i] += two._array[j+3,k]*1.3
    assert np.allclose(out._array, inp._array)
    # Blind test on all other cases
    for subscripts in others:
        inp.iadd_expand_two_to_four(subscripts, two, factor=1.7, begin0=3, end0=7, end1=4)


def test_four_index_iadd_expand_three_to_four():
    # test in detail for '1-3-1-2'
    lf = DenseLinalgFactory(4)
    inp = lf.create_four_index()
    three = lf.create_three_index()
    inp.randomize()
    out = inp.copy()
    three.randomize()
    inp.iadd_expand_three_to_four('1-3-1-2', three, factor=1.3)
    for i in range(4):
        for j in range(4):
            for k in range(4):
                out._array[i,j,i,k] += three._array[i,j,k]*1.3
    # Blind test on all other cases
    others = ['0-2-0-1', '1-2-1-2', '0-2-0-2', '0-3-0-2']
    for subscripts in others:
        inp.iadd_expand_three_to_four(subscripts, three, factor=1.7)


def test_four_index_assign_four_index_transform():
    for method in 'tensordot', 'einsum':
        lf = DenseLinalgFactory(5)
        a = lf.create_four_index()
        e0 = lf.create_expansion()
        e1 = lf.create_expansion()
        e2 = lf.create_expansion()
        e3 = lf.create_expansion()
        a.randomize()
        a.symmetrize()
        e0.randomize()
        e1.randomize()
        e2.randomize()
        e3.randomize()
        b = a.new()
        b.assign_four_index_transform(a, e0, method=method)
        assert np.allclose(b._array, b._array.transpose(1,0,3,2)), method
        assert np.allclose(b._array, b._array.transpose(2,3,0,1)), method
        assert np.allclose(b._array, b._array.transpose(1,2,3,0)), method
        assert np.isclose(b._array[0,1,2,3],
                          np.einsum('abcd,a,b,c,d->...', a._array,
                                    e0.coeffs[:,0], e0.coeffs[:,1],
                                    e0.coeffs[:,2], e0.coeffs[:,3])), method
        b.assign_four_index_transform(a, e0, e1, e2, e3, method=method)
        assert not np.allclose(b._array, b._array.transpose(1,0,3,2)), method
        assert not np.allclose(b._array, b._array.transpose(2,3,0,1)), method
        assert not np.allclose(b._array, b._array.transpose(1,2,3,0)), method
        assert np.isclose(b._array[0,1,2,3],
                          np.einsum('abcd,a,b,c,d->...', a._array,
                                    e0.coeffs[:,0], e1.coeffs[:,1],
                                    e2.coeffs[:,2], e3.coeffs[:,3])), method


def test_four_index_assign_four_index_transform_consistency():
    lf = DenseLinalgFactory(5)
    a = lf.create_four_index()
    e0 = lf.create_expansion()
    e1 = lf.create_expansion()
    e2 = lf.create_expansion()
    e3 = lf.create_expansion()
    a.randomize()
    e0.randomize()
    e1.randomize()
    e2.randomize()
    e3.randomize()
    b = a.new()
    b.assign_four_index_transform(a, e0, e1, e2, e3, method='tensordot')
    c = a.new()
    c.assign_four_index_transform(a, e0, e1, e2, e3, method='einsum')
    assert np.allclose(b._array, c._array)


#
# Tests on water (not really unit tests. oh well...)
#


def get_water_sto3g_hf(lf=None):
    if lf is None:
        lf = DenseLinalgFactory(7)
    fn = context.get_fn('test/water_sto3g_hf_g03.log')
    coeffs = np.array([
        9.94099882E-01, 2.67799213E-02, 3.46630004E-03, -1.54676269E-15,
        2.45105601E-03, -6.08393842E-03, -6.08393693E-03, -2.32889095E-01,
        8.31788042E-01, 1.03349385E-01, 9.97532839E-17, 7.30794097E-02,
        1.60223990E-01, 1.60223948E-01, 1.65502862E-08, -9.03020258E-08,
        -3.46565859E-01, -2.28559667E-16, 4.90116062E-01, 4.41542336E-01,
        -4.41542341E-01, 1.00235366E-01, -5.23423149E-01, 6.48259144E-01,
        -5.78009326E-16, 4.58390414E-01, 2.69085788E-01, 2.69085849E-01,
        8.92936017E-17, -1.75482465E-16, 2.47517845E-16, 1.00000000E+00,
        5.97439610E-16, -3.70474007E-17, -2.27323914E-17, -1.35631600E-01,
        9.08581133E-01, 5.83295647E-01, -4.37819173E-16, 4.12453695E-01,
        -8.07337352E-01, -8.07337875E-01, 5.67656309E-08, -4.29452066E-07,
        5.82525068E-01, -6.76605679E-17, -8.23811720E-01, 8.42614916E-01,
        -8.42614243E-01
    ]).reshape(7,7).T
    epsilons = np.array([
        -2.02333942E+01, -1.26583942E+00, -6.29365088E-01, -4.41724988E-01,
        -3.87671783E-01, 6.03082408E-01, 7.66134805E-01
    ])
    occ_model = AufbauOccModel(5)
    exp_alpha = lf.create_expansion()
    exp_alpha.coeffs[:] = coeffs
    exp_alpha.energies[:] = epsilons
    occ_model.assign(exp_alpha)
    assert (exp_alpha.occupations == np.array([1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0])).all()
    # convert the cache dictionary to a real cache object
    cache = Cache()
    data = load_operators_g09(fn, lf)
    for key, value in data.iteritems():
        cache.dump(key, value)
    return lf, cache, exp_alpha


def test_fock_matrix_eigen():
    lf, cache, exp_alpha = get_water_sto3g_hf()
    nbasis = cache['olp'].nbasis

    hartree = lf.create_two_index(nbasis)
    exchange = lf.create_two_index(nbasis)
    dm_alpha = exp_alpha.to_dm()
    cache['er'].contract_two_to_two('abcd,bd->ac', dm_alpha, hartree)
    cache['er'].contract_two_to_two('abcd,cb->ad', dm_alpha, exchange)

    # Construct the Fock operator
    fock = lf.create_two_index(nbasis)
    fock.iadd(cache['kin'], 1)
    fock.iadd(cache['na'], -1)
    fock.iadd(hartree, 2)
    fock.iadd(exchange, -1)

    # Check for convergence
    error = exp_alpha.error_eigen(fock, cache['olp'])
    assert error > 0
    assert error < 1e-4

    # Check self-consistency of the orbital energies
    old_energies = exp_alpha.energies.copy()
    exp_alpha.from_fock(fock, cache['olp'])
    assert abs(exp_alpha.energies - old_energies).max() < 1e-4


def test_kinetic_energy_water_sto3g():
    lf, cache, exp_alpha = get_water_sto3g_hf()
    dm = exp_alpha.to_dm()
    dm.iscale(2)
    ekin = cache['kin'].contract_two('ab,ab', dm)
    assert abs(ekin - 74.60736832935) < 1e-4


def test_ortho_water_sto3g():
    lf, cache, exp_alpha = get_water_sto3g_hf()
    for i0 in xrange(7):
        orb0 = exp_alpha.coeffs[:,i0]
        for i1 in xrange(i0+1):
            orb1 = exp_alpha.coeffs[:,i1]
            check = cache['olp'].inner(orb0, orb1)
            assert abs(check - (i0==i1)) < 1e-4


def test_potential_energy_water_sto3g_hf():
    lf, cache, exp_alpha = get_water_sto3g_hf()
    dm = exp_alpha.to_dm()
    dm.iscale(2)
    #epot = -nuclear_attraction.contract_two('ab,ab', dm)
    epot = -cache['na'].contract_two('ab,ab', dm)
    assert abs(epot - (-197.1170963957)) < 2e-3


def test_electron_electron_water_sto3g_hf():
    lf, cache, exp_alpha = get_water_sto3g_hf()
    hartree = lf.create_two_index()
    exchange = lf.create_two_index()
    dm = exp_alpha.to_dm()
    cache['er'].contract_two_to_two('abcd,bd->ac', dm, hartree)
    cache['er'].contract_two_to_two('abcd,cb->ad', dm, exchange)
    eee = 2*hartree.contract_two('ab,ab', dm) \
          - exchange.contract_two('ab,ab', dm)
    assert abs(eee - 38.29686853319) < 1e-4


def test_hartree_fock_water():
    lf, cache, exp_alpha0 = get_water_sto3g_hf()

    # Neutral water IOData
    nalpha = 5

    # Construct the hamiltonian core guess
    hamcore = lf.create_two_index()
    hamcore.iadd(cache['kin'], 1)
    hamcore.iadd(cache['na'], -1)
    exp_alpha1 = lf.create_expansion()
    exp_alpha1.from_fock(hamcore, cache['olp'])
    exp_alpha1.occupations[:nalpha] = 1.0
    assert (exp_alpha1.energies != 0.0).all()

    # The SCF loop
    hartree = lf.create_two_index()
    exchange = lf.create_two_index()
    fock = lf.create_two_index()
    for i in xrange(1000):
        # Derive the density matrix
        dm_alpha = exp_alpha1.to_dm()
        # Construct the Fock operator
        fock.clear()
        fock.iadd(hamcore, 1)
        cache['er'].contract_two_to_two('abcd,bd->ac', dm_alpha, hartree)
        cache['er'].contract_two_to_two('abcd,cb->ad', dm_alpha, exchange)
        fock.iadd(hartree, 2)
        fock.iadd(exchange, -1)
        # Check for convergence
        error = exp_alpha1.error_eigen(fock, cache['olp'])
        if error < 1e-10:
            break
        # Derive the expansion and the density matrix from the fock operator
        exp_alpha1.from_fock(fock, cache['olp'])

    assert abs(exp_alpha0.energies - exp_alpha1.energies).max() < 1e-4

    # Check the hartree-fock energy
    dm_alpha = exp_alpha1.to_dm()
    hf1 = sum([
        -2*hartree.contract_two('ab,ab', dm_alpha),
        +1*exchange.contract_two('ab,ab', dm_alpha),
    ]) + exp_alpha1.energies[:nalpha].sum()*2
    hf2 = sum([
        2*cache['kin'].contract_two('ab,ab', dm_alpha),
        -2*cache['na'].contract_two('ab,ab', dm_alpha),
        +2*hartree.contract_two('ab,ab', dm_alpha),
        -exchange.contract_two('ab,ab', dm_alpha),
    ])
    enn = 9.2535672047 # nucleus-nucleus interaction
    assert abs(hf1 + enn - (-74.9592923284)) < 1e-4
    assert abs(hf2 + enn - (-74.9592923284)) < 1e-4
