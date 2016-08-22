#!/usr/bin/env python

from horton import *


# Define Occupation model, expansion coefficients and overlap
# -----------------------------------------------------------
lf = DenseLinalgFactory(6)
occ_model = AufbauOccModel(3)
modelham = Hubbard(pbc=True)
orb = lf.create_expansion(6)
olp = modelham.compute_overlap(lf)


# One and two-body interaction terms
# ----------------------------------

# t-param, t = -1
hopping = modelham.compute_kinetic(lf, -1)
# U-param, U = 2
onsite = modelham.compute_er(lf, 2)


# Perform initial guess
# ---------------------
guess_core_hamiltonian(olp, hopping, orb)
terms = [
    RTwoIndexTerm(hopping, 'kin'),
    RDirectTerm(onsite, 'hartree'),
    RExchangeTerm(onsite, 'x_hf'),
]
ham = REffHam(terms)


# Do a Hartree-Fock calculation
# -----------------------------
scf_solver = PlainSCFSolver()
scf_solver(ham, lf, olp, occ_model, orb)


# Assign results to variables for regression testing
# --------------------------------------------------
result_energy = ham.cache["energy"]
result_orb_energies = orb.energies
result_kin = ham.cache["energy_kin"]
result_hartree = ham.cache["energy_hartree"]
result_ex = ham.cache["energy_x_hf"]
# --------------------------------------------------
