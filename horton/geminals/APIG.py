#!/usr/bin/env python2

from __future__ import division, print_function
import numpy as np
from itertools import permutations
from horton.geminals.Geminal import Geminal

class GeminalAPIG(Geminal):
    """
    Base geminal class, serving as a template for specific AP*G implementations.
    Contains information about the geminal wavefunction.

    Attributes
    ----------
    flavor : str or None
        The 'flavor' of AP*G implementation contained in the instance of the
        object.
    npairs : int
        The number of electron pairs in the geminal.
    orb : 2-index np.ndarray, dtype=np.float64
        The orbital coefficient matrix.
    coeffs : 2-index np.ndarray, dtype=np.float64
        The geminal coefficient matrix.
    proj_space : list
        A list of SlaterDets making up the projection space {\Phi}.
    overlap_cache : dict
        Cached computation results for <\Phi|\Psi> and similar.

    Methods
    -------
    clear_overlap_cache() :
        Delete the stored Slater determinant overlaps from the memory cache.
    optimize_coeffs(solver) :
        Update the coeffs matrix by minimizing {<\Phi|H|\Psi> = E<\Phi|\Psi>}
        (the projected Schrodinger equation method) according to the passed
        solver.
    optimize_orbs(optimizer) :
        Optimize the orbitals of the geminal \Psi according to the passed
        optimizer.
    compute_permanent() :
        Evaluate the permanent of the coeffs matrix.
    compute_energy() :
        Evaluate the energy of the geminal \Psi according to ham.  This can be
        simply passed to HORTON v2 as it is now.
    compute_overlap(phi) :
        Evaluate the Slater determinant overlap <phi|\Psi>.  Retrieve it from memory if it
        has previously been computed, otherwise do the computation and index
        the result in the overlap_cache.  If the overlap to be returned will not
        be accessed again, delete it from overlap_cache.
    generate_nonlinear_equations(derivative=False) :
        Compute the set of non-linear equations that can be passed to the
        projected Schrodinger equation solver.  This requires:
        * Computation of the ordered set of non-zero overlaps required for the
          equations.
        * The computation of the geminal \Psi's energy according to ham (this
          can be passed to HORTON 2.0 as it is now).
        * Optionally, the evaluation of the derivatives of the non-linear
          equations with respect to both {c} and energy.
    generate_single(phi, p, q) :
        Return the single excitation of a Slater determinant \Phi from orbital p
        to q.
    generate_double(phi, p, q, r, s) :
        Return the double excitation of a Slater determinant \Phi from orbitals
        p and q to orbitals s and r.
    """

    def __init__(self, npairs, orb, guess, proj_space=None):
        """
        Parameters
        ---------
        npairs : int
            The number of electron pairs in the geminal.
        orb : 2-index np.ndarray, dtype=np.float64
            The orbital coefficient matrix.
        guess : 2-index np.ndarray, dtype=np.float64
            The initial guess for the geminal coefficients.
        proj_space : ProjectionSpace
            A container for SlaterDets making up the projection space {\Phi} and
            cached computation results for <\Phi|\Psi> and similar.

        Raises
        ------
        NotImplementedError
        """

        self.flavor = None

        self.npairs = npairs
        self.orb = orb
        self.coeffs = guess

        if proj_space is None:
            self.proj_space = []
            for i in range(self.npairs):
                self.proj_space.append(2**i)
        else:
            self.proj_space = proj_space

        self.overlap_cache = {}

        #raise NotImplementedError

    def clear_overlap_cache(self):
        """
        Clears any stored Slater determinant overlaps from the overlap cache.
        """

        del self.overlap_cache
        self.overlap_cache = {}

    def optimize_coeffs(self, solver):
        """
        Update the coeffs matrix by minimizing {<\Phi|H|\Psi> = E<\Phi|\Psi>}
        (the projected Schrodinger equation method) according to the passed
        solver.

        Parameters
        ----------
        solver : Solver
            The callable object containing the algorithm for solving the
            projected Schrodinger equation and converging the coefficient
            matrix.
        """

        solver(self)

    def optimize_orbs(self, optimizer):
        """
        Optimize the orbitals of the geminal \Psi according to the passed
        optimizer.

        Parameters
        ----------
        optimizer : Optimizer
            The callable object containing the algorithm for optimizing the
            orbitals of the geminal \Psi.
        """

        optimizer(self)

    def compute_permanent(self, matrix=None):
        """
        Evaluate the permanent of the coeffs matrix.  The base Geminal class
        implements the *most* naive permanent evaluator using itertools
        permutations.

        Parameters
        ----------
        matrix : 2-D np.ndarray, optional
            If specified, the permanent of this matrix will be calculated
            instead of that of self.coeffs.

        Returns
        -------
        permanent :
            The permanent of `self.coeffs' or `matrix'.
        """

        if matrix is None:
            mat = self.coeffs
        else:
            mat = matrix

        assert (mat.shape[0] == mat.shape[1])

        for perm in permutations(range(mat.shape[0])):
            j, term = 0, 1
            for i in perm:
                term *= test[i][j]
                j += 1
            permanent += term

        return permanent

    def compute_energy(self, ham):
        """
        Evaluate the energy of the geminal \Psi according to ham.  This can be
        simply passed to HORTON v2 as it is now.

        Parameters
        ----------
        ham : 4-index np.ndarray, dtype=np.float64
            The Hamiltonian matrix for the geminal wavefunction.

        Returns
        -------
        Energy : np.float64
            The energy of the geminal \Psi according to ham.

        Raises
        ------
        NotImplementedError
        """

        raise NotImplementedError

    def compute_overlap(self, phi):
        """
        Evaluate the overlap <slater_det|\Psi>.  Retrieve it from memory if it
        has previously been computed, otherwise do the computation and index
        the result in the overlap_cache.  If the overlap to be returned will not
        be accessed again, delete it from overlap_cache.

        Parameters
        ----------
        phi : object
            Some sort of unique description of the Slater determinant whose
            overlap with \Psi we wish to compute.

        Returns
        -------
        overlap : object
            The result of the overlap computation for the Slater determinant phi.
        cache_flag : bool
            True if the result was pulled from cache, and False if it had to be
            computed.
        """

        if phi in self.overlap_cache:
            cache_flag = True
            print("O Cached")
            overlap = self.overlap_cache[phi]
        else:
            cache_flag = False
            print("X Computing")
            overlap  = "Fake Computation"
            self.overlap_cache[phi] = overlap

        return overlap, cache_flag

    @staticmethod
    def generate_single(phi, p, q):
        """
        Return the single excitation of a Slater determinant \Phi from orbital p
        to orbital q.

        Parameters
        ----------
        phi : int
            The Slater determinant to be excited.
        p : int
            The index of the orbital from which the excitation occurs.
        q : int
            The index of the orbital to which the excitation occurs.

        Returns
        -------
        excitation : object
            The singly-excited Slater determinant.
        """

        pmask = ~(1 << p)
        qmask = 1 << q

        excitation = phi
        excitation &= pmask
        excitation |= qmask

        return excitation

    @staticmethod
    def generate_double(phi, p, q, s, r):
        """
        Return the double excitation of a Slater determinant \Phi from orbitals
        p and q to orbitals s and r.

        Parameters
        ----------
        phi : int
            The Slater determinant to be excited.
        p : int
            The index of the orbital from which the excitation occurs.
        q : int
            The index of the orbital from which the excitation occurs.
        s : int
            The index of the orbital to which the excitation occurs.
        r : int
            The index of the orbital to which the excitation occurs.

        Returns
        -------
        excitation : object
            The doubly-excited Slater determinant.

        Raises
        ------
        AssertionError :
            Indices p and q cannot be equal, and indices s and r cannot be equal.
        """
        assert ((p != q) and (s != r))

        pmask = ~(1 << p)
        qmask = ~(1 << q)
        smask = 1 << s
        rmask = 1 << r

        excitation = phi
        excitation &= pmask
        excitation &= qmask
        excitation |= smask
        excitation |= rmask

        return excitation

    def generate_nonlinear_equations(self, ham, derivative=False):

        """
        Compute the set of non-linear equations that can be passed to the
        projected Schrodinger equation solver.  This requires:
        * Computation of the ordered set of non-zero overlaps required for the
          equations.
        * The computation of the geminal \Psi's energy according to ham (this
          can be passed to/from HORTON 2.0 as it is now).
        * Optionally, the evaluation of the derivatives of the non-linear
          equations with respect to both {c} and energy.

        Parameters
        ----------
        ham : 4-index np.ndarray, dtype=np.float64
            The Hamiltonian matrix for the geminal wavefunction.

        Returns
        -------
        equations : list of [functions | lambdas]
            A list of functions or lambdas corresponding to the non-linear
            equations to be solved.
        deriv_wrt_c : list(s) of [functions | lambdas] (if derivative=True)
            A list of 'equations' (as the above variable), where each
            corresponds to the derivative of the equation w.r.t. a
            coefficient.
        deriv_wrt_energy : list of [functions | lambdas] (if derivative=True)
            A list of 'equations' (as the above variable), where each
            corresponds to the derivative of the equation w.r.t. energy.

        Raises
        ------
        NotImplementedError
        """

        raise NotImplementedError

# -*- -*-
