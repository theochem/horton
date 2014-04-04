#!/usr/bin/env python

import os, numpy as np, sys
from glob import glob
from horton import *
from horton.meanfield.test.common import check_cubic_wrapper

#log.set_level(log.high)
log.set_level(log.silent)

debug = True

def main(fns_fchk):
    np.set_printoptions(suppress=True, linewidth=100)
    print '                            Case             Horton                G09        H-G  H converged'
    print '----------------------------------------------------------------------------------------------'
    for fn_fchk in fns_fchk:
        # Get stuff from g09
        mol = Molecule.from_file(fn_fchk)
        g09_energy = mol.energy

        # compute Gaussian integrals
        olp = mol.obasis.compute_overlap(mol.lf)
        kin = mol.obasis.compute_kinetic(mol.lf)
        na = mol.obasis.compute_nuclear_attraction(mol.coordinates, mol.pseudo_numbers, mol.lf)
        er = mol.obasis.compute_electron_repulsion(mol.lf)

        # make a list of the expansion objects
        exps = [mol.exp_alpha]
        if hasattr(mol, 'exp_beta'):
            exps.append(mol.exp_beta)

        # keep the g09 dms
        if debug:
            dms_g09 = [exp.to_dm() for exp in exps]

        # construct an initial guess
        guess_core_hamiltonian(olp, kin, na, *exps)

        # define the effective hamiltonian
        if len(exps) == 1:
            terms = [
                ROneBodyTerm(kin, 'kin'),
                RDirectTerm(er, 'hartree'),
                RExchangeTerm(er, 'x_hf'),
                ROneBodyTerm(na, 'ne'),
            ]
            ham = REffHam(terms)
        else:
            terms = [
                UOneBodyTerm(kin, 'kin'),
                UDirectTerm(er, 'hartree'),
                UExchangeTerm(er, 'x_hf'),
                UOneBodyTerm(na, 'ne'),
            ]
            ham = UEffHam(terms)

        # construct initial density matrices
        dms = [exp.to_dm() for exp in exps]

        # configure orbital occupations
        occ_model = AufbauOccModel(*[exp.occupations.sum() for exp in exps])

        # Converge the SCF
        scf_solver = ODASCFSolver(1e-8, 1024)
        niter = scf_solver(ham, mol.lf, olp, occ_model, *dms)

        # Analyze results
        horton_energy = ham.cache['energy']
        error = horton_energy - g09_energy
        if len(exps) == 1:
            prefix = 'r'
        else:
            prefix = 'u'
        print '%s %30s  % 15.10e  % 15.10e  %+9.4f  %s' % (prefix, fn_fchk, horton_energy, g09_energy, error, niter)

        if debug:
            check_cubic_wrapper(ham, dms_g09, dms, do_plot=True)


if __name__ == '__main__':
    args = sys.argv[1:]
    if len(args) == 0:
        args = sorted(glob('*/gaussian.fchk'))
    main(args)
