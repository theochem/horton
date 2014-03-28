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


from horton import *


def test_project_identical():
    sys = System.from_file(context.get_fn('test/water_sto3g_hf_g03.fchk'))
    exp = sys.lf.create_expansion()
    project_orbitals_mgs_low(sys.obasis, sys.obasis, sys.wfn.exp_alpha, exp)
    assert (exp.energies == 0.0).all()
    assert (exp.occupations == sys.wfn.exp_alpha.occupations).all()
    assert abs(exp.coeffs[:,:-2] - sys.wfn.exp_alpha.coeffs[:,:-2]).max() < 1e-9
    assert (exp.coeffs[:,-2:] == 0.0).all()


def test_project_larger():
    # Load STO3G system and keep essential results
    sys = System.from_file(context.get_fn('test/water_sto3g_hf_g03.fchk'))
    obasis0 = sys.obasis
    wfn0 = sys.wfn
    exp0 = wfn0.exp_alpha

    # Upgrade the basis to 3-21G and project
    sys.update_obasis('3-21G')
    obasis1 = sys.obasis
    setup_mean_field_wfn(sys, restricted=True)
    project_orbitals_mgs(obasis1, sys.wfn, wfn0, obasis0)
    exp1 = sys.wfn.exp_alpha
    assert (exp1.energies == 0.0).all()
    assert exp0.occupations.sum() == exp1.occupations.sum()
    assert (exp1.coeffs[:,5:] == 0.0).all()

    # Check the normalization of the projected orbitals
    olp = sys.get_overlap()
    for i0 in xrange(5):
        for i1 in xrange(i0+1):
            dot = olp.dot(sys.wfn.exp_alpha.coeffs[:,i0], sys.wfn.exp_alpha.coeffs[:,i1])
            if i0 == i1:
                assert abs(dot-1) < 1e-5
            else:
                assert abs(dot) < 1e-5

    # Setup HF hamiltonian and compute energy
    ham = Hamiltonian(sys, [HartreeFockExchange(sys.lf, sys.wfn,
                            sys.get_electron_repulsion())])
    energy1 = ham.compute()

    # Optimize wfn
    converge_scf_oda(ham)
    energy2 = sys.extra['energy']
    assert energy2 < energy1 # the energy should decrease after scf convergence

    # Construct a core initial guess
    guess_hamiltonian_core(sys)
    ham.clear()
    energy3 = ham.compute()
    assert energy3 > energy1 # the projected guess should be better than the core guess


def test_project_smaller():
    # Load 3-21G system and keep essential results
    sys = System.from_file(context.get_fn('test/li_h_3-21G_hf_g09.fchk'))
    obasis0 = sys.obasis
    wfn0 = sys.wfn

    # Downgrade the basis to sto-3g and project
    sys.update_obasis('sto-3g')
    obasis1 = sys.obasis
    setup_mean_field_wfn(sys, restricted=False)
    project_orbitals_mgs(obasis1, sys.wfn, wfn0, obasis0)
    wfn1 = sys.wfn
    assert (wfn1.exp_alpha.energies == 0.0).all()
    assert (wfn1.exp_beta.energies == 0.0).all()
    assert wfn1.exp_alpha.occupations.sum() == 2
    assert wfn1.exp_beta.occupations.sum() == 1
    assert (wfn1.exp_alpha.coeffs[:,2:] == 0.0).all()
    assert (wfn1.exp_beta.coeffs[:,1:] == 0.0).all()

    # Check the normalization of the projected orbitals
    olp = sys.get_overlap()
    for exp, nocc in (wfn1.exp_alpha, 2), (wfn1.exp_beta, 1):
        for i0 in xrange(nocc):
            for i1 in xrange(i0+1):
                dot = olp.dot(exp.coeffs[:,i0], exp.coeffs[:,i1])
                if i0 == i1:
                    assert abs(dot-1) < 1e-5
                else:
                    assert abs(dot) < 1e-5

    # Setup HF hamiltonian and compute energy
    ham = Hamiltonian(sys, [HartreeFockExchange(sys.lf, sys.wfn,
                            sys.get_electron_repulsion())])
    energy1 = ham.compute()

    # Optimize wfn
    converge_scf_oda(ham)
    energy2 = sys.extra['energy']
    assert energy2 < energy1 # the energy should decrease after scf convergence

    # Construct a core initial guess
    guess_hamiltonian_core(sys)
    ham.clear()
    energy3 = ham.compute()
    assert energy3 > energy2 # the core guess should be worse than the converged


def test_inplace():
    # Create initial system
    sys = System.from_file(context.get_fn('test/water.xyz'), obasis='sto-3g')
    setup_mean_field_wfn(sys, restricted=True)
    guess_hamiltonian_core(sys)

    old_obasis = sys.obasis
    old_wfn = sys.wfn

    # Change geometry
    sys.coordinates[1,2] += 0.1
    sys.update_obasis('sto-3g')

    # Force new wfn object
    sys._wfn = None
    setup_mean_field_wfn(sys, restricted=True)

    # Project from one to other:
    project_orbitals_mgs(sys.obasis, sys.wfn, old_wfn, old_obasis)
    exp0 = sys.wfn.exp_alpha
    assert (exp0.energies == 0.0).all()
    assert (exp0.occupations == old_wfn.exp_alpha.occupations).all()
    assert abs(exp0.coeffs[:,:5] - old_wfn.exp_alpha.coeffs[:,:5]).max() > 1e-3 # something should change
    assert (exp0.coeffs[:,5:] == 0.0).all()

    # Project in-place:
    sys._wfn = old_wfn
    project_orbitals_mgs(sys.obasis, sys.wfn, old_wfn, old_obasis)
    exp1 = sys.wfn.exp_alpha
    assert (exp1.energies == 0.0).all()
    assert (exp0.occupations == exp1.occupations).all()
    assert (exp0.coeffs == exp1.coeffs).all()
