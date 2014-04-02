#!/usr/bin/env python

import os, numpy as np, sys
from glob import glob
from horton import *
from horton.meanfield.test.common import check_cubic_os_wrapper, check_cubic_cs_wrapper

#log.set_level(log.high)
log.set_level(log.silent)

debug = False

def main(fns_fchk):
    np.set_printoptions(suppress=True, linewidth=100)
    print '                            Case             Horton                G09        H-G  H converged'
    print '----------------------------------------------------------------------------------------------'
    for fn_fchk in fns_fchk:
        # Get stuff from g09
        mol = Molecule.from_file(fn_fchk)
        if debug:
            dma0 = mol.wfn.dm_alpha.copy()
            if isinstance(mol.wfn, UnrestrictedWFN):
                dmb0 = mol.wfn.dm_beta.copy()
        g09_energy = mol.energy
        wfn = mol.wfn

        # compute Gaussian integrals
        olp = mol.obasis.compute_overlap(mol.lf)
        kin = mol.obasis.compute_kinetic(mol.lf)
        na = mol.obasis.compute_nuclear_attraction(mol.coordinates, mol.pseudo_numbers, mol.lf)
        er = mol.obasis.compute_electron_repulsion(mol.lf)

        guess_core_hamiltonian(wfn, olp, kin, na)
        if isinstance(wfn, RestrictedWFN):
            terms = [
                RestrictedOneBodyTerm(kin, 'kin'),
                RestrictedDirectTerm(er, 'hartree'),
                RestrictedExchangeTerm(er, 'x_hf'),
                RestrictedOneBodyTerm(na, 'ne'),
            ]
            ham = RestrictedEffectiveHamiltonian(terms)
        else:
            terms = [
                UnrestrictedOneBodyTerm(kin, 'kin'),
                UnrestrictedDirectTerm(er, 'hartree'),
                UnrestrictedExchangeTerm(er, 'x_hf'),
                UnrestrictedOneBodyTerm(na, 'ne'),
            ]
            ham = UnrestrictedEffectiveHamiltonian(terms)
        converged = converge_scf_oda(ham, wfn, mol.lf, olp, maxiter=1024, threshold=1e-8, debug=False)
        horton_energy = ham.cache['energy']
        error = horton_energy - g09_energy
        if isinstance(mol.wfn, RestrictedWFN):
            prefix = 'r'
        else:
            prefix = 'u'
        print '%s %30s  % 15.10e  % 15.10e  %+9.4f  %s' % (prefix, fn_fchk, horton_energy, g09_energy, error, converged)

        if debug:
            dma1 = mol.wfn.dm_alpha.copy()
            if isinstance(mol.wfn, UnrestrictedWFN):
                dmb1 = mol.wfn.dm_beta.copy()
                check_cubic_os_wrapper(ham, wfn, dma0, dmb0, dma1, dmb1, do_plot=True)
            else:
                check_cubic_cs_wrapper(ham, wfn, dma0, dma1, do_plot=True)


if __name__ == '__main__':
    args = sys.argv[1:]
    if len(args) == 0:
        args = sorted(glob('*/gaussian.fchk'))
    main(args)
