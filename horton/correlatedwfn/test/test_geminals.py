from horton import *
import numpy as np

def test_ap1rog_cs():
    fn_xyz = context.get_fn('test/h2.xyz')
    mol = Molecule.from_file(fn_xyz)
    obasis = get_gobasis(mol.coordinates, mol.numbers, '6-31G')
    lf = DenseLinalgFactory(obasis.nbasis)
    occ_model = AufbauOccModel(1)
    exp_alpha = lf.create_expansion(obasis.nbasis)
    olp = obasis.compute_overlap(lf)
    kin = obasis.compute_kinetic(lf)
    na = obasis.compute_nuclear_attraction(mol.coordinates, mol.pseudo_numbers, lf)
    guess_core_hamiltonian(olp, kin, na, exp_alpha)

    er = obasis.compute_electron_repulsion(lf)
    external = {'nn': compute_nucnuc(mol.coordinates, mol.pseudo_numbers)}
    terms = [
        RTwoIndexTerm(kin, 'kin'),
        RDirectTerm(er, 'hartree'),
        RExchangeTerm(er, 'x_hf'),
        RTwoIndexTerm(na, 'ne'),
    ]
    ham = REffHam(terms, external)
    scf_solver = PlainSCFSolver()
    scf_solver(ham, lf, olp, occ_model, exp_alpha)

    one = lf.create_two_index(obasis.nbasis)
    one.iadd(kin)
    one.iadd(na)

    # Do AP1roG optimization:
    geminal_solver = RAp1rog(lf, occ_model)
    guess = np.array([-0.08, -0.05, -0.03])
    energy, g = geminal_solver(one, er, external['nn'], exp_alpha, olp, False, **{'guess': {'geminal': guess}})
    assert (abs(energy - -1.143420629378) < 1e-6)

def test_ap1rog_cs_scf():
    fn_xyz = context.get_fn('test/h2.xyz')
    mol = Molecule.from_file(fn_xyz)
    obasis = get_gobasis(mol.coordinates, mol.numbers, '6-31G')
    lf = DenseLinalgFactory(obasis.nbasis)
    occ_model = AufbauOccModel(1)
    exp_alpha = lf.create_expansion(obasis.nbasis)
    olp = obasis.compute_overlap(lf)
    kin = obasis.compute_kinetic(lf)
    na = obasis.compute_nuclear_attraction(mol.coordinates, mol.pseudo_numbers, lf)
    guess_core_hamiltonian(olp, kin, na, exp_alpha)

    er = obasis.compute_electron_repulsion(lf)
    external = {'nn': compute_nucnuc(mol.coordinates, mol.pseudo_numbers)}
    terms = [
        RTwoIndexTerm(kin, 'kin'),
        RDirectTerm(er, 'hartree'),
        RExchangeTerm(er, 'x_hf'),
        RTwoIndexTerm(na, 'ne'),
    ]
    ham = REffHam(terms, external)
    scf_solver = PlainSCFSolver()
    scf_solver(ham, lf, olp, occ_model, exp_alpha)

    one = lf.create_two_index(obasis.nbasis)
    one.iadd(kin)
    one.iadd(na)

    # Do AP1roG optimization:
    geminal_solver = RAp1rog(lf, occ_model)
    guess = np.array([-0.1, -0.05, -0.02])
    energy, g, l = geminal_solver(one, er, external['nn'], exp_alpha, olp, True, **{'checkpoint': -1, 'guess': {'geminal': guess}})
    assert (abs(energy - -1.151686291339) < 1e-6)
