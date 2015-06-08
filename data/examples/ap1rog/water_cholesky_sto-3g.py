#! /usr/bin/python

from horton import *
import numpy as np

###############################################################################
## Set up molecule, define basis set ##########################################
###############################################################################
fn_xyz = context.get_fn('test/water.xyz')
mol = IOData.from_file(fn_xyz)
obasis = get_gobasis(mol.coordinates, mol.numbers, 'sto-3g')
###############################################################################
## Define Occupation model, expansion coefficients and overlap ################
###############################################################################
lf = CholeskyLinalgFactory(obasis.nbasis)
occ_model = AufbauOccModel(5)
orb = lf.create_expansion(obasis.nbasis)
olp = obasis.compute_overlap(lf)
###############################################################################
## Construct Hamiltonian ######################################################
###############################################################################
kin = obasis.compute_kinetic(lf)
na = obasis.compute_nuclear_attraction(mol.coordinates, mol.pseudo_numbers, lf)
er = obasis.compute_electron_repulsion(lf)
external = {'nn': compute_nucnuc(mol.coordinates, mol.pseudo_numbers)}
terms = [
    RTwoIndexTerm(kin, 'kin'),
    RDirectTerm(er, 'hartree'),
    RExchangeTerm(er, 'x_hf'),
    RTwoIndexTerm(na, 'ne'),
]
ham = REffHam(terms, external)
###############################################################################
## Perform initial guess ######################################################
###############################################################################
guess_core_hamiltonian(olp, kin, na, orb)
###############################################################################
## Do a Hartree-Fock calculation ##############################################
###############################################################################
scf_solver = PlainSCFSolver(1e-6)
scf_solver(ham, lf, olp, occ_model, orb)
###############################################################################
## Combine one-electron integrals to single Hamiltonian #######################
###############################################################################
one = kin.copy()
one.iadd(na)

###############################################################################
## Do OO-AP1roG optimization ##################################################
###############################################################################
ap1rog = RAp1rog(lf, occ_model)
energy, c, l = ap1rog(one, er, external['nn'], orb, olp, True, **{
    'indextrans': 'tensordot',
    'warning': False,
    'checkpoint': 1,
    'levelshift': 1e-8,
    'absolute': False,
    'givensrot': np.array([[]]),
    'swapa': np.array([[]]),
    'sort': True,
    'guess': {'type': 'random', 'factor': -0.1, 'geminal': None, 'lagrange': None},
    'solver': {'wfn': 'krylov', 'lagrange': 'krylov'},
    'maxiter': {'wfniter': 200, 'orbiter': 100},
    'dumpci': {'amplitudestofile': False, 'amplitudesfilename': './ap1rog_amplitudes.dat'},
    'thresh': {'wfn':  1e-12, 'energy': 1e-8, 'gradientnorm': 1e-4, 'gradientmax': 5e-5},
    'printoptions': {'geminal': True, 'ci': 0.01, 'excitationlevel': 1},
    'stepsearch': {'method': 'trust-region', 'alpha': 1.0, 'c1': 0.0001, 'minalpha': 1e-6, 'maxiterouter': 10, 'maxiterinner': 500, 'maxeta': 0.75, 'mineta': 0.25, 'upscale': 2.0, 'downscale': 0.25, 'trustradius': 0.75, 'maxtrustradius': 0.75, 'threshold': 1e-8, 'optimizer': 'ddl'},
    'orbitaloptimizer': 'variational'
}
)
