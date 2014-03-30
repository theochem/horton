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
'''A wrapper for all SCF algorithms'''


__all__ = ['SCFWrapper']


from horton.meanfield.scf import converge_scf
from horton.meanfield.scf_oda import converge_scf_oda
from horton.meanfield.scf_cdiis import converge_scf_cdiis
from horton.meanfield.scf_ediis import converge_scf_ediis
from horton.meanfield.scf_ediis2 import converge_scf_ediis2
from horton.meanfield.convergence import convergence_error_eigen, convergence_error_commutator


class SCFWrapper(object):
    '''A callable that contains pre-configure SCF options'''

    # The differen SCF routines in horton
    available_methods = {
        'plain': converge_scf,
        'oda': converge_scf_oda,
        'cdiis': converge_scf_cdiis,
        'ediis': converge_scf_ediis,
        'ediis2': converge_scf_ediis2,
    }

    # The matching convergence error functions
    error_measures = {
        'plain': convergence_error_eigen,
        'oda': convergence_error_eigen,
        'cdiis': convergence_error_commutator,
        'ediis': convergence_error_commutator,
        'ediis2': convergence_error_commutator,
    }

    def __init__(self, method, **kwargs):
        '''
           **Arguments:**

           method
                The SCF method. Select from: %s

           **Optional arguments:** depends on the SCF method. Check the
           documentation of the selected method for available options.
        ''' % (', '.join(self.available_methods))
        if method not in self.available_methods:
            raise ValueError('Unknown SCF method: %s' % method)
        self.method = method
        self.kwargs = kwargs

    def __call__(self, ham, wfn, lf, overlap):
        '''Converge the SCF for the given Hamiltonian

           ham
                A Hamiltonian instance

           **Returns:** the number of iterations
        '''
        return self.available_methods[self.method](ham, wfn, lf, overlap, **self.kwargs)

    def convergence_error(self, ham, wfn, lf, overlap):
        return self.error_measures[self.method](ham, wfn, lf, overlap)
