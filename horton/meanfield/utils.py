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
"""Utility functions"""


__all__ = [
    'check_dm', 'get_level_shift', 'get_spin', 'get_homo_lumo',
    'compute_commutator',
]


def check_dm(dm, overlap, lf, eps=1e-4, occ_max=1.0):
    '''Check if the density matrix has eigenvalues in the proper range.

       **Arguments:**

       dm
            The density matrix

       overlap
            The overlap matrix

       lf
            A LinalgFactory instance.

       **Optional arguments:**

       eps
            The threshold on the eigenvalue inequalities.

       occ_max
            The maximum occupation.

       A ValueError is raised when the density matrix has illegal eigenvalues.
    '''
    # construct natural orbitals
    exp = lf.create_expansion()
    exp.derive_naturals(dm, overlap)
    if exp.occupations.min() < -eps:
        raise ValueError('The density matrix has eigenvalues considerably smaller than zero. error=%e' % (exp.occupations.min()))
    if exp.occupations.max() > occ_max+eps:
        raise ValueError('The density matrix has eigenvalues considerably larger than one. error=%e' % (exp.occupations.max()-1))


def get_level_shift(dm, overlap):
    '''Construct a level shift operator.

       **Arguments:**

       dm
            A density matrix.

       overlap
            The overlap matrix

       **Returns:** The level-shift operator.
    '''
    level_shift = overlap.copy()
    level_shift.idot(dm)
    level_shift.idot(overlap)
    return level_shift


def get_spin(exp_alpha, exp_beta, overlap):
    '''Returns the expectation values of the projected and squared spin

       **Arguments:**

       exp_alpha, exp_beta
            The alpha and beta orbitals.

       overlap
            The overlap matrix

       **Returns:** sz, ssq
    '''
    nalpha = exp_alpha.occupations.sum()
    nbeta = exp_beta.occupations.sum()
    sz = (nalpha - nbeta)/2
    correction = 0.0
    for ialpha in xrange(exp_alpha.nfn):
        if exp_alpha.occupations[ialpha] == 0.0:
            continue
        for ibeta in xrange(exp_beta.nfn):
            if exp_beta.occupations[ibeta] == 0.0:
                continue
            correction += overlap.inner(exp_alpha.coeffs[:,ialpha],
                                        exp_beta.coeffs[:,ibeta])**2

    ssq = sz*(sz+1) + nbeta - correction
    print sz, ssq
    return sz, ssq


def get_homo_lumo(*exps):
    '''Return the HOMO and LUMO energy for the given expansion

       **Arguments:**

       exp1, exp2, ...
            DensityExpansion objects

       **Returns:** homo_energy, lumo_energy. (The second is None when all
       orbitals are occupied.)
    '''
    homo_energy = max(exp.homo_energy for exp in exps)
    lumo_energies = [exp.lumo_energy for exp in exps]
    lumo_energies = [lumo_energy for lumo_energy in lumo_energies if lumo_energy is not None]
    if len(lumo_energies) == 0:
        lumo_energy = None
    else:
        lumo_energy = min(lumo_energies)
    return homo_energy, lumo_energy


def compute_commutator(dm, fock, overlap, work, output):
    '''Compute the dm-fock commutator, including an overlap matrix

       **Arguments:** (all TwoIndex objects)

       dm
            A density matrix

       fock
            A fock matrix

       overlap
            An overlap matrix

       work
            A temporary matrix

       output
            The output matrix in which the commutator, S.D.F-F.D.S, is stored.
    '''
    # construct sdf
    work.assign(overlap)
    work.idot(dm)
    work.idot(fock)
    output.assign(work)
    # construct fds and subtract
    work.assign(fock)
    work.idot(dm)
    work.idot(overlap)
    output.iadd(work, factor=-1)
