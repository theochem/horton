#!/usr/bin/env python

from horton import *
import h5py as h5
import os

log.set_level(log.silent)


def store_wfn(fn_h5, mixing, name_case, exp):
    with h5.File(fn_h5) as f:
        name_mixing = '%08.5f' % (-np.log10(mixing))
        grp = f.require_group(name_mixing)
        grp = grp.require_group(name_case)

        # clear the group if anything was present
        for key in list(grp.keys()):
            del grp[key]
        for key in list(grp.attrs.keys()):
            del grp.attrs[key]
        exp.to_hdf5(grp)
        # The following is needed to create object of the right type when
        # reading from the checkpoint:
        grp.attrs['class'] = exp.__class__.__name__


def get_random_occupations(nbasis, nep):
    result = np.zeros(nbasis)
    # this is not uniformely random, but it is good enough.
    for iep in range(int(np.round(nep))):
        total = 1.0
        while total > 0:
            if total < 0.01:
                fraction = total
                total = 0.0
            else:
                fraction = np.random.uniform(0, total)
                total -= fraction
            index = np.random.randint(nbasis)
            result[index] += fraction
            if result[index] > 1:
                total += result[index] - 1
                result[index] = 1.0
    return result

def main():
    try:
        os.remove("guesses.h5")
    except OSError:
        pass

    fn_name = context.get_fn('test/2h-azirine.xyz')

    mol = Molecule.from_file(fn_name)
    obasis = get_gobasis(mol.coordinates, mol.numbers, '3-21G')
    lf = DenseLinalgFactory(obasis.nbasis)

    # Compute Gaussian integrals
    olp = obasis.compute_overlap(lf)
    kin = obasis.compute_kinetic(lf)
    na = obasis.compute_nuclear_attraction(mol.coordinates, mol.pseudo_numbers, lf)
    er = obasis.compute_electron_repulsion(lf)

    # Create alpha orbitals
    exp_alpha = lf.create_expansion()

    # Initial guess
    guess_core_hamiltonian(olp, kin, na, exp_alpha)

    # Construct the restricted HF effective Hamiltonian
    external = {'nn': compute_nucnuc(mol.coordinates, mol.pseudo_numbers)}
    terms = [
        ROneBodyTerm(kin, 'kin'),
        RDirectTerm(er, 'hartree'),
        RExchangeTerm(er, 'x_hf'),
        ROneBodyTerm(na, 'ne'),
    ]
    ham = REffHam(terms, external)

    # Decide how to occupy the orbitals (5 alpha electrons)
    occ_model = AufbauOccModel(5)

    # Converge WFN with plain SCF
    scf_solver = PlainSCFSolver(1e-6)
    scf_solver(ham, lf, olp, occ_model, exp_alpha)

    # generate randomized wavefunctions:
    # - arbitrary unitary transformation
    # - arbitrary (fractional) occupation numbers (with proper sum)
    nbasis = obasis.nbasis
    random_exps = []
    for irandom in range(nrandom):
        # random symmetric matrix
        tmp1 = np.random.normal(0, 1, (nbasis, nbasis))
        tmp1 = tmp1 + tmp1.T
        # the random unitary matrix
        utrans = np.linalg.eigh(tmp1)[1]
        # apply transformation
        coeffs = np.dot(exp_alpha.coeffs, utrans)
        # random occupation numbers
        occupations = get_random_occupations(nbasis, exp_alpha.occupations.sum())
        # create a expansion object
        exp_alpha_temp = lf.create_expansion()
        # assign the random orbitals
        exp_alpha_temp.coeffs[:] = coeffs
        exp_alpha_temp.occupations[:] = occupations
        # store the expansion in the h5 file and in the list
        store_wfn('guesses.h5', 1.0, 'case_%03i' % irandom, exp_alpha_temp)
        random_exps.append(exp_alpha_temp)

    # interpolate between solution and random wfns
    for mixing in mixings[1:]: # do not consider mixing==1.0
        for irandom in range(nrandom):
            # create a new wfn object.
            # construct the mixed density matrix
            dm_mixed = lf.create_one_body()
            dm_mixed.iadd(random_exps[irandom].to_dm(), mixing)
            dm_mixed.iadd(ham.cache['dm_alpha'], 1-mixing)
            # turn it into a set of orbitals
            exp_alpha_temp = lf.create_expansion()
            exp_alpha_temp.derive_naturals(dm_mixed, olp)
            # store the wfn in the h5 file

            store_wfn('guesses.h5', mixing, 'case_%03i' % irandom, exp_alpha_temp)


if __name__ == '__main__':
    mixings = np.array([1, 1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8])
    nrandom = 20
    main()
