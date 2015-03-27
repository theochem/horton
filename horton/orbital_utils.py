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
'''Utility functions for orbital modifications'''


from horton.log import timer
from horton.matrix import TwoIndex, Expansion


__all__ = ['rotate_orbitals',
           'compute_unitary_matrix',
           'transform_integrals',
           'read_orbitals',
          ]



def rotate_orbitals(*args):
    '''Rotate AO/MO (or MO/MO) coefficient matrix such that C` = C*U

       **Arguments:**

       args
            Rotation matrices (TwoIndex instance) and AO/MO coefficients
            (Expansion instance). rots and exps are ordered such that rot1
            corresponds to exp1, etc., i.e., rot1, rot2,..., exp1, exp2,...
    '''
    exps = []
    rotmat = []
    for arg in args:
        if isinstance(arg, TwoIndex):
            rotmat.append(arg)
        elif isinstance(arg, Expansion):
            exps.append(arg)
        else:
            raise TypeError('argument of unsupported type: %s' % arg)
    for i in xrange(len(exps)):
        exps[i].assign_dot(exps[i], rotmat[i])


def compute_unitary_matrix(kappa):
    '''Determine a unitary matrix from a skew-symmetric matrix K as
       U = exp(-K) by approximating U = 1 - K + 1/2 K^2 + O(3)

       **Arguments:**

       kappa
            A skew-symmetric matrix (TwoIndex instance)
    '''
    out = kappa.new()
    out.assign_diagonal(1.0)

    #
    # Approximate unitary matrix
    # U = exp(-K) by U = 1 - K + 1/2 K^2 + O(3)
    #
    out.iadd(kappa, -1.0)
    out.iadd_dot(kappa, kappa, 0.5)
    #
    # orthogonalization because approximate U matrix might not be unitary/orthogonal:
    #
    out.iortho()
    return out


@timer.with_section('Index Trans')
def transform_integrals(one, two, indextrans='tensordot', *exps):
    '''Update MO integrals. Returns list of transformed 1- and 2-electron
       integrals according to a list of expansion coefficients.

       **Arguments:**

       one
           One-electron integrals in the AO basis. A TwoIndex instance.

       two
           Two-electron integrals in the AO basis.

       **Optional arguments:**

       indextrans
           Choice of 4-index transformation. Default 'tensordot'.

       args
           The expansion coefficients.

    '''
    exp = []
    two_mo = []
    one_mo = []
    if not all([isinstance(i, Expansion) for i in exps]):
        raise TypeError('argument of unsupported type: %s' %i)
    for arg in exps:
        exp.append(arg)

    #
    # Loop over all possible AO/MO coefficients. We distinguish between
    #   * restricted orbitals: only one expansion instance (alpha=beta), returns
    #                          one set of 1- and 2-body integrals
    #   * unrestricted orbitals: two expansion instances (alpha, beta), returns
    #                            two sets of 1-body, and three sets of 2-body
    #                            integrals (one_alpha,alpha and one_beta,beta)
    #                            and (<alpha alpha|alpha alpha>, )
    #                            (<alpha beta|alpha beta>, <beta beta|beta beta>)
    #                            Note that order of alpha and beta is determined
    #                            by the order of the expansion instances
    #
    for i in xrange(len(exps)):
        for j in xrange(i, len(exps)):
            #
            # Transform 2-electron part
            # Note that integrals are stored using physics' convention <12|12>
            #
            out4ind = two.new()
            out4ind.assign_four_index_transform(two, exp[i], exp[j], exp[i], exp[j], indextrans)
            two_mo.append(out4ind)

        #
        # Transform 1-electron part
        #
        out2ind = one.new()
        out2ind.assign_two_index_transform(one, exps[i])
        one_mo.append(out2ind)
    return one_mo, two_mo


def read_orbitals(exp, olp, orbfile="./orb.hdf5", olpfile="./olp.hdf5"):
    '''Read orbitals from hdf5-file. Function can be used if orbitals from
       different geometries are used as initial guess.

       **Arguments:**

       exp
            The wfn expansion coefficients (Expansion instance)

       olp
            The AO overlap matrix (TwoIndex instance)

       orbfile
            The filename of the orbitals (str)

       olpfile
            The filename of the overlap matrix (str)
    '''
    import h5py
    forb = h5py.File(orbfile, "r")
    folp = h5py.File(olpfile, "r")
    newexp = exp.from_hdf5(forb)
    newolp = olp.from_hdf5(folp)
    #
    # Calculate Sr^{1/2}
    #
    newolp12 = newolp.sqrt()
    #
    # Get square root of current AO overlap matrix
    #
    olp12 = olp.sqrt()
    #
    # Calculate S^{-1/2}
    #
    olpinv = olp12.inverse()
    #
    # Calculate Sr^{1/2}*C
    #
    newolp12.idot(newexp)
    #
    # Get new AO/MO coefficient matrix from S^{-1/2}*Sr^{1/2}*C that
    # satisfies C^T*S*C=1
    #
    olpinv.idot(newolp12)
    exp.assign(olpinv)

    forb.close()
    folp.close()
