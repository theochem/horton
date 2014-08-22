#!/usr/bin/env python

from horton import *
import h5py as h5
import os
from scoop import futures
import scipy as scp

log.set_level(log.silent)

class MultipleFailedConvergence(Exception):
    pass

class Method(object):
    def __init__(self):
        raise NotImplementedError
    def solve(self):
        raise NotImplementedError

# class DirectOpt(Method):
#     def __init__(self, lf, ham):
#         fock = lf.create_one_body()
#         ham.compute_fock(fock, None)
#         sys.wfn.clear() #FIXME: update for new structure after guess code is ported.
#         sys.wfn.update_exp(fock, sys.get_overlap())
#         args, N = initialGuess.setup_guess(sys, ifCheat=True, isFrac=True)
#         cons = initialGuess.setup_cons(sys, args, N, restricted=True)
#         self.lg,self.x0 = initialGuess.setup_lg(sys, ham, cons, args, "None", "None", isFrac=True, restricted=True)
#     def solve(self):
#         x_star, niter = newton_rewrite.solve(self.lg, self.x0)
#         return niter

class SCFOpt(Method):
    def __init__(self, ham, lf, olp, occ_model, dm_alpha):
        self.ham = ham
        self.scf_solver = EDIIS2SCFSolver(1e-8)
        self.lf = lf
        self.olp = olp
        self.occ_model = occ_model
        self.dm = dm_alpha
        
    def solve(self):
        return self.scf_solver(self.ham, self.lf, self.olp, self.occ_model, self.dm)

def load_exp(fn_h5, mixing, name_case, lf):
    with h5.File(fn_h5) as f:
        name_mixing = '%08.5f' % (-np.log10(mixing))
        group_mixing = f[name_mixing]
        group_case = group_mixing[name_case]
        # cheating is nice ...
        class_name = group_case.attrs.get('class')
        cls = __import__('horton', fromlist=[class_name]).__dict__[class_name]
        return cls.from_hdf5(group_case, lf)

def run(irandom, mixing, method):
    nfail = 0
    mix_iters = []

    fn_name = context.get_fn('test/2h-azirine.xyz')

    mol = Molecule.from_file(fn_name)
    obasis = get_gobasis(mol.coordinates, mol.numbers, '3-21G')
    lf = DenseLinalgFactory(obasis.nbasis)
    
    
    olp = obasis.compute_overlap(lf)
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
    
    # Decide how to occupy the orbitals (5 alpha electrons)
    occ_model = AufbauOccModel(5)
    
    for i in mixing[::-1]:
        exp_alpha = load_exp('guesses.h5', i, 'case_%03i' % irandom, lf)
        print "mixing: " + str(i)
        
        dm_alpha = exp_alpha.to_dm()
        ham.reset(dm_alpha)
        
        energy0 = ham.compute()
        solver = method(ham, lf, olp, occ_model, dm_alpha)
        
        try:
            niter = solver.solve()
            mix_iters.append(niter)
        except (scp.optimize.nonlin.NoConvergence,NoSCFConvergence):
            nfail += 1
            if nfail > 2:
                padding = ["inf"]*(len(mixing) - len(mix_iters))
                return mix_iters + padding
    
        print 'run %3i: %8.5f %3i %12.6f %12.6f' % (
            irandom, -np.log10(i), niter, energy0, ham.compute()
            )
    return mix_iters

def main():
    try:
        os.remove("scf_results.txt")
    except OSError:
        pass
#     
#     res =  run(5, mixings, SCFOpt)
#     print [str(mix) +" " + str(iter) for mix,iter in zip(reversed(mixings[-len(res):]), res)]
    
    fseq = [futures.submit(run, i, mixings, SCFOpt) for i in xrange(nrandom)]
          
    not_done = ["dummy"]
    while not_done:
        done, not_done = futures.wait(fseq, None, "FIRST_COMPLETED")
              
        for i in done:
            with open('scf_results.txt', 'a') as f:
                line = [str(mix) +" "+ str(iter) for mix,iter in zip(reversed(mixings[-len(i.result()):]), i.result())]
                print >> f, '\n'.join(line)
            fseq.remove(i)
    
if __name__ == '__main__':
    mixings = np.array([1, 1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8])
    nrandom = 20
    main()
