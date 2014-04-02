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
"""Occupation number models"""


import numpy as np

from horton.log import log
from horton.exceptions import ElectronCountError
from horton.quadprog import find_1d_root
from horton.constants import boltzmann


__all__ = [
    'AufbauOccModel', 'AufbauSpinOccModel', 'FermiOccModel',
]


class AufbauOccModel(object):
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

    def assign(self, *exps):
        '''apply occupation numbers to the expansion objects

           **Arguments:**

           exp_alpha, exp_beta, ...
                Expansion objects
        '''
        if len(exps) != len(self.noccs):
            raise TypeError('Expected %i expansion objects, got %i.' % (len(self.nocc), len(exps)))
        for exp, nocc in zip(exps, self.noccs):
            if exp.nfn < nocc:
                raise ElectronCountError('The number of orbitals must not be lower than the number of alpha or beta electrons.')
            # It is assumed that the orbitals are sorted from low to high energy.
            if nocc == int(nocc):
                exp.occupations[:nocc] = 1.0
                exp.occupations[nocc:] = 0.0
            else:
                exp.occupations[:int(np.floor(nocc))] = 1.0
                exp.occupations[int(np.floor(nocc))] = nocc - np.floor(nocc)
                exp.occupations[int(np.ceil(nocc)):] = 0.0

    def log(self):
        log('Occupation model: %s' % self)
        log.deflist([('nalpha', self.nalpha), ('nbeta', self.nbeta)])

    def check_dms(self, overlap, *dms, **kwargs):
        eps = kwargs.pop('eps', 1e-4)
        if len(kwargs) > 0:
            raise TypeError('Unexpected keyword arguments: %s' % kwargs.keys())
        for dm, nocc in zip(dms, self.noccs):
            assert abs(overlap.expectation_value(dm) - nocc) < eps


class AufbauSpinOccModel(object):
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

    def assign(self, exp_alpha, exp_beta):
        '''apply occupation numbers to the expansion objects

           **Arguments:**

           exp_alpha, exp_beta, ...
                Expansion objects
        '''
        nel = self.nel
        ialpha = 0
        ibeta = 0
        while nel > 0:
            if exp_alpha.energies[ialpha] <= exp_beta.energies[ibeta]:
                exp_alpha.occupations[ialpha] = min(1.0, nel)
                ialpha += 1
            else:
                exp_beta.occupations[ibeta] = min(1.0, nel)
                ibeta += 1
            nel -= 1

    def log(self):
        log('Occupation model: %s' % self)
        log.deflist([('nel', self.nel)])

    def check_dms(self, overlap, *dms, **kwargs):
        eps = kwargs.pop('eps', 1e-4)
        if len(kwargs) > 0:
            raise TypeError('Unexpected keyword arguments: %s' % kwargs.keys())
        assert abs(sum(overlap.expectation_value(dm) for dm in dms) - self.nel) < eps


class FermiBase(object):
    '''Base class for Fermi smearing occupation models'''
    def __init__(self, temperature=300, eps=1e-8):
        '''
           **Optional arguments:**

           temperature
                Controls the width of the distribution (derivative)

           eps
                The error on the sum of the occupation number when searching for
                the right Fermi level.
        '''
        if temperature <= 0:
            raise ValueError('The temperature must be strictly positive')
        if eps <= 0:
            raise ValueError('The root-finder threshold (eps) must be strictly positive.')
        self.temperature = float(temperature)
        self.eps = eps


class FermiOccModel(FermiBase, AufbauOccModel):
    '''Fermi smearing electron occupation model'''
    def __init__(self, *noccs, **kwargs):
        '''
           **Arguments:**

           nalpha, nbeta, ...
                The number of electrons in each channel.

           **Optional keyword arguments:**

           temperature
                Controls the width of the distribution (derivative)

           eps
                The error on the sum of the occupation number when searching for
                the right Fermi level.
        '''
        temperature = kwargs.pop('temperature', 300)
        eps = kwargs.pop('eps', 1e-8)
        if len(kwargs) > 0:
            raise TypeError('Unknown keyword arguments: %s' % kwargs.keys())
        FermiBase.__init__(self, temperature, eps)
        AufbauOccModel.__init__(self, *noccs)

    def assign(self, *exps):
        '''apply occupation numbers to the expansion objects

           **Arguments:**

           exp_alpha, exp_beta, ...
                Expansion objects
        '''
        beta = 1.0/self.temperature/boltzmann
        for exp, nocc in zip(exps, self.noccs):
            def get_occ(mu):
                occ = np.zeros(exp.nfn)
                mask = exp.energies < mu
                e = np.exp(beta*(exp.energies[mask] - mu))
                occ[mask] = 1.0/(e + 1.0)
                mask = ~mask
                e = np.exp(-beta*(exp.energies[mask] - mu))
                occ[mask] = e/(1.0 + e)
                return occ

            def error(mu):
                return nocc - get_occ(mu).sum()

            mu0 = exp.energies[exp.nfn/2]
            error0 = error(mu0)
            delta = 0.1*(1 - 2*(error0 < 0))
            for i in xrange(100):
                mu1 = mu0 + delta
                error1 = error(mu1)
                if error1 == 0 or ((error0 > 0) ^ (error1 > 0)):
                    break
                delta *= 2

            if error1 == 0:
                exp.occupations[:] = get_occ(mu1)
            else:
                mu, error = find_1d_root(error, (mu0, error0), (mu1, error1), eps=self.eps)
                exp.occupations[:] = get_occ(mu)

    def log(self):
        log('Occupation model: %s' % self)
        log.deflist([
            ('nalpha', self.nalpha),
            ('nbeta', self.nbeta),
            ('temperature', self.temperature),
            ('eps', self.eps),
        ])
