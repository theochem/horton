# -*- coding: utf-8 -*-
# Horton is a Density Functional Theory program.
# Copyright (C) 2011-2012 Toon Verstraelen <Toon.Verstraelen@UGent.be>
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


import numpy as np
from horton import *


def test_hamiltonian_init():
    coordinates = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
    numbers = np.array([1, 1])
    sys = System(coordinates, numbers, obasis='STO-3G')
    sys.init_wfn(0, 1)

    # test if no terms gives a ValueError
    try:
        ham = Hamiltonian(sys, [])
        assert False
    except ValueError:
        pass

    # test if terms are added automagically
    ham = Hamiltonian(sys, [HartreeFock()])
    assert sum(isinstance(term, KineticEnergy) for term in ham.terms) == 1
    assert sum(isinstance(term, ExternalPotential) for term in ham.terms) == 1
    assert sum(isinstance(term, Hartree) for term in ham.terms) == 1

    # test if the necessary operators are constructed in the system object
    assert 'kin' in sys.operators
    assert 'na' in sys.operators
    assert 'er' in sys.operators
    assert 'olp' in sys.operators

    # check attribute of HartreeFock term
    assert ham.terms[0].fraction_exchange == 1.0


def test_energy_hydrogen():
    fn_fchk = context.get_fn('test/h_sto3g.fchk')
    sys = System.from_file(fn_fchk)
    ham = Hamiltonian(sys, [HartreeFock()])
    ham.compute_energy()
    assert abs(sys.props['energy'] - -4.665818503844346E-01) < 1e-8
