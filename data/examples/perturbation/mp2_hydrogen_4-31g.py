from horton import *
###############################################################################
## Set up molecule, define basis set ##########################################
###############################################################################
mol = Molecule.from_file('hydrogen.xyz')
obasis = get_gobasis(mol.coordinates, mol.numbers, '4-31G')
###############################################################################
## Define Occupation model, expansion coefficients and overlap ################
###############################################################################
#lf = DenseLinalgFactory(obasis.nbasis)
lf = CholeskyLinalgFactory(obasis.nbasis)
occ_model = AufbauOccModel(1)
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
## Get Hartree-Fock energy ####################################################
###############################################################################
ehf = ham.compute()
###############################################################################
## Combine one-electron integrals to single Hamiltonian #######################
###############################################################################
one = kin.copy()
one.iadd(na)

###############################################################################
## Do RMP2 calculation ########################################################
###############################################################################
mp2 = RMP2(lf, occ_model)
emp2, tmp2 = mp2(one, er, orb, **{'eref': ehf})
