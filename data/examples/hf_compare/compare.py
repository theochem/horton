#!/usr/bin/env python

from glob import glob
import numpy as np
import os
import sys

from horton import *  # pylint: disable=wildcard-import,unused-wildcard-import
from horton.meanfield.test.common import check_cubic_wrapper


# log.set_level(log.high)
log.set_level(log.silent)

debug = True


def main(fns_fchk):
    """Main program.

    Parameters
    ----------
    fns_fchk : str
        A list of paths to formatted checkpoint files to be reproduced.
    """
    np.set_printoptions(suppress=True, linewidth=100)
    print('                            Case             HORTON                G09        H-G  H converged')
    print('----------------------------------------------------------------------------------------------')
    for fn_fchk in fns_fchk:
        # Get stuff from g09.
        mol = IOData.from_file(fn_fchk)
        g09_energy = mol.energy

        # Compute Gaussian integrals.
        olp = mol.obasis.compute_overlap()
        kin = mol.obasis.compute_kinetic()
        na = mol.obasis.compute_nuclear_attraction(mol.coordinates, mol.pseudo_numbers)
        er = mol.obasis.compute_electron_repulsion()

        # Make a list of the expansion objects.
        orbs = [mol.orb_alpha]
        if hasattr(mol, 'orb_beta'):
            orbs.append(mol.orb_beta)

        # Keep the g09 dms.
        if debug:
            dms_g09 = [orb.to_dm() for orb in orbs]

        # Construct an initial guess.
        guess_core_hamiltonian(olp, kin + na, *orbs)

        # Define the effective hamiltonian.
        if len(orbs) == 1:
            terms = [
                RTwoIndexTerm(kin, 'kin'),
                RDirectTerm(er, 'hartree'),
                RExchangeTerm(er, 'x_hf'),
                RTwoIndexTerm(na, 'ne'),
            ]
            ham = REffHam(terms)
        else:
            terms = [
                UTwoIndexTerm(kin, 'kin'),
                UDirectTerm(er, 'hartree'),
                UExchangeTerm(er, 'x_hf'),
                UTwoIndexTerm(na, 'ne'),
            ]
            ham = UEffHam(terms)

        # Construct initial density matrices.
        dms = [orb.to_dm() for orb in orbs]

        # Configure orbital occupations.
        noccs = np.round(np.array([orb.occupations.sum() for orb in orbs])).astype(int)
        occ_model = AufbauOccModel(*noccs)

        # Converge the SCF.
        scf_solver = ODASCFSolver(1e-8, 1024)
        niter = scf_solver(ham, olp, occ_model, *dms)

        # Analyze results.
        horton_energy = ham.cache['energy']
        error = horton_energy - g09_energy
        if len(orbs) == 1:
            prefix = 'r'
        else:
            prefix = 'u'
        print('%s %30s  % 15.10e  % 15.10e  %+9.4f  %s' % (
            prefix, fn_fchk, horton_energy, g09_energy, error, niter))

        if debug:
            check_cubic_wrapper(ham, dms_g09, dms, do_plot=True)


if __name__ == '__main__':
    args = sys.argv[1:]
    if len(args) == 0:
        args = sorted(glob('*/gaussian.fchk'))
    main(args)
