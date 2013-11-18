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


import numpy as np, tempfile, shutil
from contextlib import contextmanager
import subprocess, os, shlex

from horton import *


__all__ = [
    'check_script', 'check_delta',
    'get_random_cell', 'get_pentagon_moments',
    'compare_expansions', 'compare_all_expansions', 'compare_dms',
    'compare_all_dms', 'compare_one_body', 'compare_two_body',
    'compare_occ_model', 'compare_wfns', 'compare_systems',
    'compare_symmetries',
    'tmpdir',
]


# All, except underflows, is *not* fine.
np.seterr(divide='raise', over='raise', invalid='raise')


def check_script(command, workdir):
    env = dict(os.environ)
    root_dir = os.getcwd()
    env['PYTHONPATH'] = root_dir + ':' + env.get('PYTHONPATH', '')
    if os.path.isfile('setup.py') and os.path.isdir('horton') and os.path.isdir('data'):
        # This is only needed when running the tests from the source tree.
        env['HORTONDATA'] = os.path.join(root_dir, 'data')
    env['PATH'] = os.path.join(root_dir, 'scripts') + ':' + env.get('PATH', '')
    proc = subprocess.Popen(shlex.split(command), stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=workdir, env=env)
    outdata, errdata = proc.communicate()
    if proc.returncode != 0:
        print 'Standard output'
        print '+'*80
        print outdata
        print '+'*80
        print 'Standard error'
        print '+'*80
        print errdata
        print '+'*80
        assert False


def check_delta(fun, fun_deriv, x, dxs):
    """Check the difference between two function values using the analytical gradient

       Arguments:

       fun
            The function whose derivatives must be to be tested

       fun_deriv
            The implementation of the analytical derivatives

       x
            The argument for the reference point.

       dxs
            A list with small relative changes to x

       For every displacement in ``dxs``, the following computation is repeated:

       1) D1 = 'fun(x+dx) - fun(x)' is computed.
       2) D2 = '0.5 (fun_deriv(x+dx) + fun_deriv(x)) . dx' is computed.

       A threshold is set to the median of the D1 set. For each case where |D1|
       is larger than the threshold, |D1 - D2|, should be smaller than the
       threshold.
    """
    if len(dxs) < 20:
        raise ValueError('At least 20 displacements are needed for good statistics.')

    dn1s = []
    dn2s = []
    dnds = []
    f0 = fun(x)
    grad0 = fun_deriv(x)
    for dx in dxs:
        f1 = fun(x+dx)
        grad1 = fun_deriv(x+dx)
        grad = 0.5*(grad0+grad1)
        d1 = f1 - f0
        if hasattr(d1, '__iter__'):
            norm = np.linalg.norm
        else:
            norm = abs
        d2 = np.dot(grad.ravel(), dx.ravel())

        dn1s.append(norm(d1))
        dn2s.append(norm(d2))
        dnds.append(norm(d1-d2))
    dn1s = np.array(dn1s)
    dn2s = np.array(dn2s)
    dnds = np.array(dnds)

    # Get the threshold (and mask)
    threshold = np.median(dn1s)
    mask = dn1s > threshold
    # Make sure that all cases for which dn1 is above the treshold, dnd is below
    # the threshold
    if not (dnds[mask] < threshold).all():
        raise AssertionError((
            'The first order approximation on the difference is too wrong. The '
            'threshold is %.1e.\n\nDifferences:\n%s\n\nFirst order '
            'approximation to differences:\n%s\n\nAbsolute errors:\n%s')
            % (threshold,
            ' '.join('%.1e' % v for v in dn1s[mask]),
            ' '.join('%.1e' % v for v in dn2s[mask]),
            ' '.join('%.1e' % v for v in dnds[mask])
        ))


def get_random_cell(a, nvec):
    if nvec == 0:
        return Cell(None)
    while True:
        if a <= 0:
            raise ValueError('The first argument must be strictly positive.')
        rvecs = np.random.uniform(0, a, (nvec,3))
        cell = Cell(rvecs)
        if cell.volume > a**nvec*0.1:
            return cell


def get_pentagon_moments(rmat=None, lmax=4):
    if rmat is None:
        rmat = np.identity(3, float)

    cartesian_powers = get_cartesian_powers(lmax)
    ncart = cartesian_powers.shape[0]
    result = np.zeros(ncart)
    for i in xrange(6):
        alpha = 2.0*np.pi/5.0
        vec = np.array([1+np.cos(alpha), np.sin(alpha), 0])
        vec = np.dot(rmat, vec)
        for i in xrange(ncart):
            px, py, pz = cartesian_powers[i]
            result[i] = vec[0]**px + vec[1]**py + vec[2]**pz
    return result


def compare_expansions(wfn1, wfn2, spin):
    if 'exp_%s' % spin in wfn1._cache:
        assert 'exp_%s' % spin in wfn2._cache
        e1 = wfn1.get_exp(spin)
        e2 = wfn2.get_exp(spin)
        assert e1.nbasis == e2.nbasis
        assert e1.nfn == e2.nfn
        assert (e1.coeffs == e2.coeffs).all()
        assert (e1.energies == e2.energies).all()
        assert (e1.occupations == e2.occupations).all()
    else:
        assert 'exp_%s' % spin not in wfn2._cache


def compare_all_expansions(wfn1, wfn2):
    compare_expansions(wfn1, wfn2, 'alpha')
    compare_expansions(wfn1, wfn2, 'beta')


def compare_dms(wfn1, wfn2, select):
    if 'dm_%s' % select in wfn1._cache:
        assert 'dm_%s' % select in wfn2._cache
        dm1 = wfn1.get_dm(select)
        dm2 = wfn2.get_dm(select)
        assert dm1.nbasis == dm2.nbasis
        assert (dm1._array == dm2._array).all()
    else:
        assert 'dm_%s' % select not in wfn2._cache


def compare_all_dms(wfn1, wfn2):
    compare_dms(wfn1, wfn2, 'alpha')
    compare_dms(wfn1, wfn2, 'beta')
    compare_dms(wfn1, wfn2, 'full')
    compare_dms(wfn1, wfn2, 'spin')


def compare_one_body(sys1, sys2, key):
    if key in sys1.cache:
        assert key in sys2.cache
        op1 = sys1.cache[key]
        op2 = sys2.cache[key]
        if isinstance(op1, DenseOneBody):
            assert isinstance(op2, DenseOneBody)
            assert op1.nbasis == op2.nbasis
            assert (op1._array == op2._array).all()
        else:
            raise NotImplementedError
    else:
        assert key not in sys2.cache


def compare_two_body(sys1, sys2, key):
    if key in sys1.cache:
        assert key in sys2.cache
        op1 = sys1.cache[key]
        op2 = sys2.cache[key]
        if isinstance(op1, DenseTwoBody):
            assert isinstance(op2, DenseTwoBody)
            assert op1.nbasis == op2.nbasis
            assert (op1._array == op2._array).all()
        else:
            raise NotImplementedError
    else:
        assert key not in sys2.cache


def compare_occ_model(occ_model1, occ_model2):
    assert occ_model1.__class__ == occ_model2.__class__
    if occ_model1 is None:
        assert occ_model2 is None
    elif isinstance(occ_model1, AufbauOccModel):
        assert occ_model1.nalpha == occ_model2.nalpha
        assert occ_model1.nbeta == occ_model2.nbeta
    else:
        raise NotImplementedError


def compare_wfns(wfn1, wfn2):
    if isinstance(wfn1, RestrictedWFN):
        assert isinstance(wfn2, RestrictedWFN)
        assert wfn1.nbasis == wfn2.nbasis
        assert wfn1.norb == wfn2.norb
        compare_all_expansions(wfn1, wfn2)
        compare_all_dms(wfn1, wfn2)
        compare_occ_model(wfn1.occ_model, wfn2.occ_model)
    elif isinstance(wfn1, UnrestrictedWFN):
        assert isinstance(wfn2, UnrestrictedWFN)
        assert wfn1.nbasis == wfn2.nbasis
        assert wfn1.norb == wfn2.norb
        compare_all_expansions(wfn1, wfn2)
        compare_all_dms(wfn1, wfn2)
        compare_occ_model(wfn1.occ_model, wfn2.occ_model)
    elif wfn1 is None:
        assert wfn2 is None
    else:
        raise NotImplementedError


def compare_systems(sys1, sys2):
    assert (sys1.numbers == sys2.numbers).all()
    assert (sys1.coordinates == sys2.coordinates).all()
    # orbital basis
    if sys1.obasis is not None:
        assert (sys1.obasis.centers == sys2.obasis.centers).all()
        assert (sys1.obasis.shell_map == sys2.obasis.shell_map).all()
        assert (sys1.obasis.nprims == sys2.obasis.nprims).all()
        assert (sys1.obasis.shell_types == sys2.obasis.shell_types).all()
        assert (sys1.obasis.alphas == sys2.obasis.alphas).all()
        assert (sys1.obasis.con_coeffs == sys2.obasis.con_coeffs).all()
    else:
        assert sys2.obasis is None
    # wfn
    compare_wfns(sys1.wfn, sys2.wfn)
    # one-body operators
    compare_one_body(sys1, sys2, 'olp')
    compare_one_body(sys1, sys2, 'kin')
    compare_one_body(sys1, sys2, 'na')
    # two-body operators
    compare_two_body(sys1, sys2, 'er')


def compare_symmetries(s0, s1):
    assert s0.name == s1.name
    assert (s0.generators == s1.generators).all()
    assert (s0.fracs == s1.fracs).all()
    assert (s0.numbers == s1.numbers).all()
    assert s0.cell.nvec == s1.cell.nvec
    assert (s0.cell.rvecs == s1.cell.rvecs).all()
    assert (s0.labels == s1.labels).all()


@contextmanager
def tmpdir(name):
    dn = tempfile.mkdtemp(name)
    try:
        yield dn
    finally:
        shutil.rmtree(dn)
