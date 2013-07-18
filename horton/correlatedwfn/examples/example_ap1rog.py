'''
This is an example input file to run an optimization
of AP1roG with and without orbital optimization.

In the following, all input parameters are specified
using their default values.

A minimal input example is also given.
'''
from horton import *

###########################################################################################
## Set up molecule, define basis set ######################################################
###########################################################################################
fn_xyz = context.get_fn('test/h2.xyz')
mol = Molecule.from_file(fn_xyz)
obasis = get_gobasis(mol.coordinates, mol.numbers, 'cc-pvdz')
###########################################################################################
## Define Occupation model, expansion coefficients and overlap ############################
###########################################################################################
lf = DenseLinalgFactory(obasis.nbasis)
occ_model = AufbauOccModel(1)
exp_alpha = lf.create_expansion(obasis.nbasis)
olp = obasis.compute_overlap(lf)
###########################################################################################
## Construct Hamiltonian ##################################################################
###########################################################################################
kin = obasis.compute_kinetic(lf)
na = obasis.compute_nuclear_attraction(mol.coordinates, mol.pseudo_numbers, lf)
er = obasis.compute_electron_repulsion(lf)
external = {'nn': compute_nucnuc(mol.coordinates, mol.pseudo_numbers)}
terms = [
    ROneBodyTerm(kin, 'kin'),
    RDirectTerm(er, 'hartree'),
    RExchangeTerm(er, 'x_hf'),
    ROneBodyTerm(na, 'ne'),
]
ham = REffHam(terms, external)
###########################################################################################
## Perform initial guess ##################################################################
###########################################################################################
guess_core_hamiltonian(olp, kin, na, exp_alpha)
###########################################################################################
## Do a Hartree-Fockk calculation #########################################################
###########################################################################################
scf_solver = PlainSCFSolver()
scf_solver(ham, lf, olp, occ_model, exp_alpha)

###########################################################################################
## Combine to single one-electron Hamiltonian #############################################
###########################################################################################
one = lf.create_one_body(obasis.nbasis)
one.iadd(kin)
one.iadd(na)

###########################################################################################
## Optimize AP1roG; minimal input example #################################################
###########################################################################################
geminal_solver = RAp1rog(lf, occ_model)
energy, g = geminal_solver(one, er, external['nn'], exp_alpha, olp, False)

###########################################################################################
## Input with all default function parameters: ############################################
###########################################################################################
## Each parameter can be changed separately ###############################################
###########################################################################################
energy, g = geminal_solver(one, er, external['nn'], exp_alpha, olp, False,
                           **{'indextrans': 'tensordot',
                              'warning': False,
                              'swapa': [-1,-1],
                              'guess': {'guess': 'random', 'factor': -0.1},
                              'solver': {'wfn': 'krylov'},
                              'maxiter': {'wfn': 128},
                              'dumpci': {'amplitudestofile': False, 'amplitudesfilename': "./ap1rog_amplitudes.dat"},
                              'thresh': {'wfn': 1e-12},
                              'printoptions': {'geminal': True, 'ci': 0.01, 'excitationlevel': 1},
                             }
                           )
###########################################################################################
## Do OO-AP1roG optimization: #############################################################
## Minimal input (default parameters are set to sensible values): #########################
############################This is all you need!!!########################################
energy, g, l = geminal_solver(one, er, external['nn'], exp_alpha, olp, True)
###########################################################################################
# Input with all default function parameters:        ######################################
###########################################################################################
## Each parameter can be changed separately ###############################################
###########################################################################################
energy, g, l = geminal_solver(one, er, external['nn'], exp_alpha, olp, True,
                              **{'indextrans': 'tensordot',
                                 'warning': False,
                                 'checkpoint': 1,
                                 'lshift': 1e-8,
                                 'swapa': [-1,-1],
                                 'givensrotationa': [-1,-1,0],
                                 'guess': {'guess': 'random', 'factor': -0.1},
                                 'solver': {'wfn': 'krylov', 'lagrange': 'krylov'},
                                 'maxiter': {'wfn': 128, 'orbiter': 50},
                                 'dumpci': {'amplitudestofile': False, 'amplitudesfilename': "./ap1rog_amplitudes.dat"},
                                 'thresh': {'wfn': 1e-12, 'energy': 1e-8, 'gradientnorm': 1e-8, 'gradientmax': 5e-5},
                                 'printoptions': {'geminal': True, 'ci': 0.01, 'excitationlevel': 1},
                                 'stepsearch': {'method': 'trust-region', 'stepa': 1.0, 'c1': 0.0001, 'c2': 0.9,
                                                'maxstep': 0.75, 'minstep': 1e-6, 'maxiterouter': 40, 'maxeta': 0.5, 'optimizer': 'pcg',
                                                'maxtrustradius': 0.75, 'mineta': 0.0, 'upscale': 1.2, 'downscale': 0.7,
                                                'trustradius': 0.75, 'threshold': 1e-8, 'maxiterinner': 500}
                                 }
                              )
###########################################################################################
