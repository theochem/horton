# -*- coding: utf-8 -*-
# HORTON: Helpful Open-source Research TOol for N-fermion systems.
# Copyright (C) 2011-2022 The HORTON Development Team
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
"""Occupation number models"""

import numpy as np

from horton.exceptions import ElectronCountError
from horton.quadprog import find_1d_root
from horton.constants import boltzmann
from horton.log import biblio
from horton.utils import doc_inherit

__all__ = [
    'FixedOccModel', 'AufbauOccModel', 'AufbauSpinOccModel', 'FermiOccModel',
]


class OccModel(object):
    '''Base class for the occupation models'''

    def assign(self, *orbs):
        '''Assign occupation numbers to the Orbitals objects

           **Arguments:**

           orb_alpha, orb_beta, ...
                Orbitals objects
        '''
        raise NotImplementedError

    def check_dms(self, overlap, *dms, **kwargs):
        '''Test if the given density matrices contain the right number of electrons

           **Arguments:**

           overlap
                The overlap operator.

           dm1, dm2, ...
                Density matrices to be tested.

           **Optional keyword arguments:**

           eps (default=1e-4)
                The allowed deviation.
        '''
        raise NotImplementedError


class FixedOccModel(OccModel):
    def __init__(self, *occ_arrays):
        self.occ_arrays = occ_arrays

    @doc_inherit(OccModel)
    def assign(self, *orbs):
        if len(orbs) != len(self.occ_arrays):
            raise TypeError('Expected %i Orbitals objects, got %i.' % (len(self.nocc), len(orbs)))
        for orb, occ_array in zip(orbs, self.occ_arrays):
            orb.occupations[:len(occ_array)] = occ_array
            orb.occupations[len(occ_array):] = 0.0

    @doc_inherit(OccModel)
    def check_dms(self, overlap, *dms, **kwargs):
        eps = kwargs.pop('eps', 1e-4)
        if len(kwargs) > 0:
            raise TypeError('Unexpected keyword arguments: %s' % list(kwargs.keys()))
        if len(dms) != len(self.occ_arrays):
            raise TypeError('The number of density matrices is incorrect.')
        for dm, occ_array in zip(dms, self.occ_arrays):
            assert abs(np.einsum('ab,ba', overlap, dm) - occ_array.sum()) < eps


class AufbauOccModel(OccModel):
    '''The standard Aufbau occupation number model.

       This model just fills up all the lowest lying orbitals. When the total
       number of electrons in one channel is fractional, the fractional electron
       is put in the HOMO orbital.
    '''

    def __init__(self, *noccs):
        '''
           **Arguments:**

           nalpha, nbeta, ...
                The number of electrons in each channel.
        '''
        for nocc in noccs:
            if nocc < 0:
                raise ElectronCountError('Negative number of electrons is not allowed.')
        if sum(noccs) == 0:
            raise ElectronCountError('At least one electron is required.')
        self.noccs = noccs

    @doc_inherit(OccModel)
    def assign(self, *orbs):
        if len(orbs) != len(self.noccs):
            raise TypeError('Expected %i Orbitals objects, got %i.' % (len(self.nocc), len(orbs)))
        for orb, nocc in zip(orbs, self.noccs):
            if orb.nfn < nocc:
                raise ElectronCountError(
                    'The number of orbitals must not be lower than the number of alpha or beta electrons.')
            # It is assumed that the orbitals are sorted from low to high energy.
            if nocc == int(nocc):
                orb.occupations[:nocc] = 1.0
                orb.occupations[nocc:] = 0.0
            else:
                orb.occupations[:int(np.floor(nocc))] = 1.0
                orb.occupations[int(np.floor(nocc))] = nocc - np.floor(nocc)
                orb.occupations[int(np.ceil(nocc)):] = 0.0

    @doc_inherit(OccModel)
    def check_dms(self, overlap, *dms, **kwargs):
        eps = kwargs.pop('eps', 1e-4)
        if len(kwargs) > 0:
            raise TypeError('Unexpected keyword arguments: %s' % list(kwargs.keys()))
        if len(dms) != len(self.noccs):
            raise TypeError('The number of density matrices is incorrect.')
        for dm, nocc in zip(dms, self.noccs):
            assert abs(np.einsum('ab,ba', overlap, dm) - nocc) < eps


class AufbauSpinOccModel(OccModel):
    '''This Aufbau model only applies to unrestricted wavefunctions'''

    def __init__(self, nel):
        '''
           **Arguments:**

           nel
                The total number of electrons (alpha + beta)
        '''
        if nel <= 0:
            raise ElectronCountError('The number of electron must be positive.')
        self.nel = nel

    @doc_inherit(OccModel)
    def assign(self, orb_alpha, orb_beta):
        nel = self.nel
        ialpha = 0
        ibeta = 0
        while nel > 0:
            if orb_alpha.energies[ialpha] <= orb_beta.energies[ibeta]:
                orb_alpha.occupations[ialpha] = min(1.0, nel)
                ialpha += 1
            else:
                orb_beta.occupations[ibeta] = min(1.0, nel)
                ibeta += 1
            nel -= 1

    @doc_inherit(OccModel)
    def check_dms(self, overlap, *dms, **kwargs):
        eps = kwargs.pop('eps', 1e-4)
        if len(kwargs) > 0:
            raise TypeError('Unexpected keyword arguments: %s' % list(kwargs.keys()))
        assert abs(sum(np.einsum('ab,ba', overlap, dm) for dm in dms) - self.nel) < eps


class FermiOccModel(AufbauOccModel):
    '''Fermi smearing electron occupation model'''

    def __init__(self, *noccs, **kwargs):
        r'''
           **Arguments:**

           nalpha, nbeta, ...
                The number of electrons in each channel.

           **Optional keyword arguments:**

           temperature
                Controls the width of the distribution (derivative)

           eps
                The error on the sum of the occupation number when searching for
                the right Fermi level.

           For each channel, the orbital occupations are assigned with the Fermi
           distribution:

           .. math::

                n_i = \frac{1}{1 + e^{(\epsilon_i - \mu)/k_B T}}

           where, for a given set of energy levels, :math:`\{\epsilon_i\}`, the
           chemical potential, :math:`\mu`, is optimized as to satisfy the
           following constraint:

           .. math::

               \sum_i n_i = n_\text{occ}

           where :math:`n_\text{occ}` can be set per (spin) channel. This is
           only a part of the methodology presented in [rabuck1999]_.
        '''
        temperature = kwargs.pop('temperature', 300)
        eps = kwargs.pop('eps', 1e-8)
        if len(kwargs) > 0:
            raise TypeError('Unknown keyword arguments: %s' % list(kwargs.keys()))
        if temperature <= 0:
            raise ValueError('The temperature must be strictly positive')
        if eps <= 0:
            raise ValueError('The root-finder threshold (eps) must be strictly positive.')
        self.temperature = float(temperature)
        self.eps = eps
        AufbauOccModel.__init__(self, *noccs)
        biblio.cite('rabuck1999', 'the Fermi broading method to assign orbital occupations')

    @doc_inherit(OccModel)
    def assign(self, *orbs):
        beta = 1.0 / self.temperature / boltzmann
        for orb, nocc in zip(orbs, self.noccs):
            def get_occ(mu):
                occ = np.zeros(orb.nfn)
                mask = orb.energies < mu
                e = np.exp(beta * (orb.energies[mask] - mu))
                occ[mask] = 1.0 / (e + 1.0)
                mask = ~mask
                e = np.exp(-beta * (orb.energies[mask] - mu))
                occ[mask] = e / (1.0 + e)
                return occ

            def error(mu):
                return nocc - get_occ(mu).sum()

            mu0 = orb.energies[orb.nfn // 2]
            error0 = error(mu0)
            delta = 0.1 * (1 - 2 * (error0 < 0))
            for i in range(100):
                mu1 = mu0 + delta
                error1 = error(mu1)
                if error1 == 0 or ((error0 > 0) ^ (error1 > 0)):
                    break
                delta *= 2

            if error1 == 0:
                orb.occupations[:] = get_occ(mu1)
            else:
                mu, error = find_1d_root(error, (mu0, error0), (mu1, error1), eps=self.eps)
                orb.occupations[:] = get_occ(mu)
