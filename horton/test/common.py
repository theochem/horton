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


from contextlib import contextmanager
import numpy as np
import os
import shutil
import subprocess
import tempfile

from horton.cext import Cell
from horton.moments import get_cartesian_powers
from horton.matrix.dense import DenseTwoIndex, DenseFourIndex
from horton.meanfield.occ import AufbauOccModel


__all__ = [
    'in_horton_source_root',
    'check_script', 'check_script_in_tmp', 'check_delta',
    'get_random_cell', 'get_pentagon_moments',
    'compare_expansions', 'compare_all_expansions', 'compare_dms',
    'compare_all_dms', 'compare_operators', 'compare_occ_model', 'compare_exps',
    'compare_mols', 'compare_symmetries',
    'tmpdir', 'numpy_seed', 'truncated_file',
]


# All, except underflows, is *not* fine.
np.seterr(divide='raise', over='raise', invalid='raise')


def in_horton_source_root():
    '''Test if the current directory is the HORTON source tree root'''
    # Check for some files and directories that must be present (for functions
    # that use this check).
    if not os.path.isfile('setup.py'):
        return False
    if not os.path.isdir('horton'):
        return False
    if not os.path.isdir('data'):
        return False
    if not os.path.isdir('scripts'):
        return False
    if not os.path.isfile('README'):
        return False
    with open('README') as f:
        if f.next() != 'HORTON: *H*elpful *O*pen-source *R*esearch *TO*ol for *N*-fermion systems.\n':
            return False
    return True


def check_script(command, workdir):
    '''Change to workdir and try to run the given command.

       This can be used to test whether a script runs properly, with exit code
       0. When the example generates an output file, then use
       :py:function:`horton.test.common.check_script_in_tmp` instead.

       **Arguments:**

       command
            The command to be executed as a single string.

       workdir
            The work directory where the command needs to be executed.
    '''
    env = dict(os.environ)
    root_dir = os.getcwd()
    env['PYTHONPATH'] = root_dir + ':' + env.get('PYTHONPATH', '')
    if in_horton_source_root():
        # This is only needed when running the tests from the source tree.
        env['HORTONDATA'] = os.path.join(root_dir, 'data')
    env['PATH'] = os.path.join(root_dir, 'scripts') + ':' + env.get('PATH', '')
    try:
        proc = subprocess.Popen(command, stdin=subprocess.PIPE,
                                stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                cwd=workdir, env=env, shell=True)
        outdata, errdata = proc.communicate()
    except OSError:
        raise AssertionError('Executable not found.')
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


def check_script_in_tmp(command, required, expected):
    '''Test a script in a tmp dir

       **Arguments:**

       command
            The command to be executed in a tmp dir.

       required
            The required files to be copied to the tmp dir.

       expected
            A list of files expected to be present in the tmp dir after
            execution.
    '''
    with tmpdir('check_scrip_in_tmp-%s' % os.path.basename(command)) as dn:
        # copy files into tmp
        for fn in required:
            shutil.copy(fn, os.path.join(dn, os.path.basename(fn)))
        # run the script
        check_script(command, dn)
        # check the output files
        for fn in expected:
            assert os.path.exists(os.path.join(dn, fn)), 'Missing output %s' % fn


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
    assert len(x.shape) == 1
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
        d2 = np.dot(grad, dx)

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
        for j in xrange(ncart):
            px, py, pz = cartesian_powers[j]
            result[j] += vec[0]**px * vec[1]**py * vec[2]**pz
    return result


def get_point_moments(coordinates, rmat=None, lmax=4):
    if rmat is None:
        rmat = np.identity(3, float)

    cartesian_powers = get_cartesian_powers(lmax)
    ncart = cartesian_powers.shape[0]
    result = np.zeros(ncart)
    for i in xrange(len(coordinates)):
        vec = np.dot(rmat, coordinates[i])
        for j in xrange(ncart):
            px, py, pz = cartesian_powers[j]
            result[j] += vec[0]**px * vec[1]**py * vec[2]**pz
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


def compare_operators(op1, op2):
    if isinstance(op1, DenseTwoIndex) or isinstance(op1, DenseFourIndex):
        assert isinstance(op2, op1.__class__)
        assert op1.nbasis == op2.nbasis
        assert (op1._array == op2._array).all()
    else:
        raise NotImplementedError


def compare_occ_model(occ_model1, occ_model2):
    assert occ_model1.__class__ == occ_model2.__class__
    if occ_model1 is None:
        assert occ_model2 is None
    elif isinstance(occ_model1, AufbauOccModel):
        assert occ_model1.nalpha == occ_model2.nalpha
        assert occ_model1.nbeta == occ_model2.nbeta
    else:
        raise NotImplementedError


def compare_exps(exp1, exp2):
    assert exp1.nbasis == exp2.nbasis
    assert exp1.nfn == exp2.nfn
    assert (exp1.coeffs == exp2.coeffs).all()
    assert (exp1.energies == exp2.energies).all()
    assert (exp1.occupations == exp2.occupations).all()


def compare_mols(mol1, mol2):
    assert (getattr(mol1, 'title') == getattr(mol2, 'title'))
    assert (mol1.numbers == mol2.numbers).all()
    assert (mol1.coordinates == mol2.coordinates).all()
    # orbital basis
    if mol1.obasis is not None:
        assert (mol1.obasis.centers == mol2.obasis.centers).all()
        assert (mol1.obasis.shell_map == mol2.obasis.shell_map).all()
        assert (mol1.obasis.nprims == mol2.obasis.nprims).all()
        assert (mol1.obasis.shell_types == mol2.obasis.shell_types).all()
        assert (mol1.obasis.alphas == mol2.obasis.alphas).all()
        assert abs(mol1.obasis.con_coeffs - mol2.obasis.con_coeffs).max() < 1e-8
    else:
        assert mol2.obasis is None
    # wfn
    for key in 'exp_alpha', 'exp_beta':
        if hasattr(mol1, key):
            assert hasattr(mol2, key)
            compare_exps(getattr(mol1, key), getattr(mol2, key))
        else:
            assert not hasattr(mol2, key)
    # operators
    for key in 'olp', 'kin', 'na', 'er', 'dm_full_mp2', 'dm_spin_mp2', \
               'dm_full_mp3', 'dm_spin_mp3', 'dm_full_ci', 'dm_spin_ci', \
               'dm_full_cc', 'dm_spin_cc', 'dm_full_scf', 'dm_spin_scf':
        if hasattr(mol1, key):
            assert hasattr(mol2, key)
            compare_operators(getattr(mol1, key), getattr(mol2, key))
        else:
            assert not hasattr(mol2, key)


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


@contextmanager
def numpy_seed(seed=1):
    """Temporarily set NumPy's random seed to a given number.

    Parameters
    ----------
    seed : int
           The seed for NumPy's random number generator.
    """
    state = np.random.get_state()
    np.random.seed(seed)
    yield None
    np.random.set_state(state)


@contextmanager
def truncated_file(name, fn_orig, nline, nadd):
    """Make a temporary truncated copy of a file.

    Parameters
    ----------
    name : str
           The name of test, used to make a unique temporary directory
    fn_orig : str
              The file to be truncated.
    nline : int
            The number of lines to retain.
    nadd : int
           The number of empty lines to add.
    """
    with tmpdir(name) as dn:
        fn_truncated = '%s/truncated_%i_%s' % (dn, nline, os.path.basename(fn_orig))
        with open(fn_orig) as f_orig, open(fn_truncated, 'w') as f_truncated:
            for counter, line in enumerate(f_orig):
                if counter >= nline:
                    break
                f_truncated.write(line)
            for _ in xrange(nadd):
                f_truncated.write('\n')
        yield fn_truncated
