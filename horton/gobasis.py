# -*- coding: utf-8 -*-
# Horton is a Density Functional Theory program.
# Copyright (C) 2011-2012 Toon Verstraelen <Toon.Verstraelen@UGent.be>
#
# This file is part of Horton.
#
# Horton is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
#
# Horton is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, see <http://www.gnu.org/licenses/>
#
#--
"""Gaussian orbital basis set module."""


import numpy as np

from horton.context import context
from horton.periodic import periodic
from horton.cext import get_shell_nbasis


__all__ = [
    'GOBasisDesc', 'str_to_shell_types', 'GOBasisFamily', 'go_basis_families',
    'GOBasisAtom', 'GOBasisContraction', 'GOBasis'
]


class GOBasisDesc(object):
    """A user specification of the basis set."""
    def __init__(self, default, element_map=None, index_map=None, pure=True):
        """
           **Arguments:**

           default
                The default basis set applied to each atom.

           **Optional arguments:**

           element_map
                A dictionary with element names or numbers as keys, and
                basis sets as values. These specs override the default basis

           index_map
                A dictionary with atomic indexes (based on the order of the
                atoms) as keys and basis sets as values

           pure
                By default pure basis functions are used. Set this to false to
                switch to Cartesian basis functions.

           In all three arguments, basis may either refer to None, a string
           representation (STO-3G), a GOBasisFamily instance or a GOBasisAtom
           instance.
        """
        self.default = default
        if element_map is None:
            self.element_map = {}
        else:
            self.element_map = element_map
        if index_map is None:
            self.index_map = {}
        else:
            self.index_map = index_map
        self.pure = pure

        # Update the element map such that it only contains numbers as keys.
        for key in self.element_map.keys():
            if isinstance(key, basestring):
                number = periodic[key].number
                self.element_map[number] = element_map[key]
                del element_map[key]

    def apply_to(self, system):
        """Construct a GOBasis object for the system"""
        shell_map = []
        nprims = []
        shell_types = []
        alphas = []
        con_coeffs = []

        def get_basis(i, n):
            """Look up the basis for a given atom"""
            basis = self.index_map.get(i)
            if basis is not None:
                return basis
            basis = self.element_map.get(n)
            if basis is not None:
                return basis
            if self.default is None:
                raise KeyError('Could not find basis for atom %i.' % i)
            else:
                return self.default

        def translate_basis(basis_x, n):
            """Translate the first argument into a GOBasisAtom instance"""
            if isinstance(basis_x, basestring):
                basis_fam = go_basis_families.get(basis_x.lower())
                if basis_fam is None:
                    raise ValueError('Unknown basis family: %s' % basis_x)
                basis_atom = basis_fam.get(n)
                if basis_atom is None:
                    raise ValueError('The basis family %s does not contain element %i.' % (basis_x, n))
                return basis_atom
            elif isinstance(basis_x, GOBasisFamily):
                basis_atom = basis_x.get(n)
                if basis_atom is None:
                    raise ValueError('The basis family %s does not contain element %n.' % (basis_x.name, n))
                return basis_atom
            elif isinstance(basis_x, GOBasisAtom):
                return basis_x
            else:
                raise ValueError('Can not interpret %s as an atomic basis function.' % basis_x)

        # Loop over the atoms and fill in all the lists
        for i in xrange(system.natom):
            n = system.numbers[i]
            basis_x = get_basis(i, n)
            basis_atom = translate_basis(basis_x, n)
            basis_atom.extend(i, shell_map, nprims, shell_types, alphas, con_coeffs, self.pure)

        # Turn the lists into arrays
        shell_map = np.array(shell_map)
        nprims = np.array(nprims)
        shell_types = np.array(shell_types)
        alphas = np.array(alphas)
        con_coeffs = np.array(con_coeffs)

        # Return the Gaussian basis object.
        return GOBasis(shell_map, nprims, shell_types, alphas, con_coeffs)


def str_to_shell_types(s):
    """Convert a string into a list of contraction types"""
    d = {'s': 0, 'p': 1, 'd': 2, 'f': 3, 'g': 4, 'h': 5, 'i': 6}
    return [d[c] for c in s.lower()]


class GOBasisFamily(object):
    def __init__(self, name, basis_atom_map=None, filename=None):
        """
           **Arguments:**

           name
                The name of the basis set family, e.g. 'STO-3G'.

           **Optional Arguments:**

           basis_atom_map
                A dictionary with atomic numbers as values and GOBasisAtom
                instances as keys.

           filename
                A file to load the basis set from when needed.

           Either one of the two optional arguments must be provided. The first
           allows the user to craft his/her own basis set family. The second
           can be used to refer to a file that contains the basis set family.
           This file is only loaded when a basis set is requested for a given
           atom, in order to avoid slow startup times. The format of the basis
           set filename is deduced from the extension.
        """
        if basis_atom_map is None and filename is None:
            raise ValueError('Either one of the two optional arguments must be provided.')
        self.name = name
        self.basis_atom_map = basis_atom_map
        self.filename = filename

    def get(self, number):
        """Return the GOBasisAtom instance for the given number."""
        if self.basis_atom_map is None:
            self.load()
        return self.basis_atom_map.get(number)

    def load(self):
        """Load the basis set from file."""
        if self.filename.endswith('.nwchem'):
            self._load_nwchem()
        else:
            raise IOError('File format not supported: %s' % self.filename)

    def _load_nwchem(self):
        '''Load the basis set family from an NWChem file.'''
        # TODO: move this to io.py
        f = open(self.filename)
        self.basis_atom_map = {}
        bc = None # The current contraction being loaded
        for line in f:
            # strip of comments and white space
            line = line[:line.find('#')].strip()
            if len(line) == 0:
                continue
            if line == 'END':
                break
            if line.startswith('BASIS'):
                continue
            words = line.split()
            if words[0].isalpha():
                # A new contraction begins, maybe even a new atom.
                n = periodic[words[0]].number
                ba = self.basis_atom_map.get(n)
                shell_types = str_to_shell_types(words[1])
                bcs = [GOBasisContraction(shell_type, [], []) for shell_type in shell_types]
                if ba is None:
                    ba = GOBasisAtom(bcs)
                    self.basis_atom_map[n] = ba
                else:
                    ba.bcs.extend(bcs)
            else:
                # An extra primitive for the current contraction(s).
                exponent = float(words[0])
                coeffs = [float(w) for w in words[1:]]
                for i, bc in enumerate(bcs):
                    bc.alphas.append(exponent)
                    bc.con_coeffs.append(coeffs[i::len(bcs)])
        f.close()
        self._to_arrays()
        self._to_segmented()

    def _to_arrays(self):
        '''Convert all contraction attributes to numpy arrays'''
        for ba in self.basis_atom_map.itervalues():
            for bc in ba.bcs:
                bc.to_arrays()

    def _to_segmented(self):
        '''Convert all contractions from generalized to segmented'''
        new_basis_atom_map = {}
        for n, ba in self.basis_atom_map.iteritems():
            new_bcs = []
            for bc in ba.bcs:
                new_bcs.extend(bc.get_segmented_bcs())
            new_ba = GOBasisAtom(new_bcs)
            new_basis_atom_map[n] = new_ba
        self.basis_atom_map = new_basis_atom_map



go_basis_families = [
    GOBasisFamily('STO-3G', filename=context.get_fn('basis/sto-3g.nwchem')),
    GOBasisFamily('3-21G', filename=context.get_fn('basis/3-21g.nwchem')),
    GOBasisFamily('6-31++G**', filename=context.get_fn('basis/6-21++g**.nwchem')),
    GOBasisFamily('ANO', filename=context.get_fn('basis/ano-rcc.nwchem')),
]
go_basis_families = dict((bf.name.lower(), bf) for bf in go_basis_families)


class GOBasisAtom(object):
    """Description of an atomic basis set with segmented contractions."""
    def __init__(self, bcs):
        """
           **Arguments:**

           bcs
                A list of GOBasisContraction instances.
        """
        self.bcs = bcs

    def extend(self, i, shell_map, nprims, shell_types, alphas, con_coeffs, pure=True):
        """Extend the lists with basis functions for this atom.

           **Arguments:**

           i
                The index of the center of this atom.

           shell_map, nprims, shell_types, alphas, con_coeffs
                These are all output arguments that must be extended with the
                basis set parameters for this atom. The meaning of each argument
                is defined in the documentation for constructor of the GOBasis
                class.

           **Optional argument:**

           pure
                By default pure basis functions are used. Set this to False to
                switch to Cartesian basis functions.
        """
        for bc in self.bcs:
            if bc.is_generalized():
                raise ValueError('Generalized contractions are not supported (yet).')
            shell_map.append(i)
            nprims.append(len(bc.alphas))
            if pure or bc.shell_type < 2:
                shell_types.append(bc.shell_type)
            else:
                shell_types.append(-bc.shell_type)
            alphas.extend(bc.alphas)
            con_coeffs.extend(bc.con_coeffs)


class GOBasisContraction(object):
    """Description of a contraction"""
    def __init__(self, shell_type, alphas, con_coeffs):
        """
           **Arguments:**

           shell_type
                The angular quantum numbers of the shell. (0=s, 1=p, 2=d, ...)

           alphas
                1D array with exponents of the primitives.

           con_coeffs
                the linear coefficients of the contraction. In the case of a
                segmented basis set, this is just a 1D array with the same size
                as alphas. In the case of a generalized contraction, this is a
                2D numpy array, where each row corresponds to a primitive and
                the columns correspond to different contractions.

           It is possible to construct this object with lists instead of arrays
           as arguments. Just call the to_arrays method once the lists are
           completed. (This is convenient when loading basis sets from a file.)
        """
        self.shell_type = shell_type
        self.alphas = alphas
        self.con_coeffs = con_coeffs

    def to_arrays(self):
        '''Convert the alphas and con_coeffs attributes to numpy arrays.'''
        self.alphas = np.array(self.alphas)
        self.con_coeffs = np.array(self.con_coeffs)

    def __length__(self):
        '''Return the length of the contraction'''
        l = len(self.alphas)
        assert l == len(self.con_coeffs)
        return l

    def is_generalized(self):
        '''Returns True if this is a generalized contraction

           This routine also checks if the con_coeffs attribute is consistent.
        '''
        if len(self.con_coeffs.shape) == 2:
            return True
        elif len(self.con_coeffs.shape) == 1:
            return False
        else:
            raise ValueError

    def get_segmented_bcs(self):
        '''Returns a list of segmented contractions'''
        if not self.is_generalized():
            raise TypeError('Conversion to segmented contractions only makes sense for generalized contractions.')
        return [
            GOBasisContraction(self.shell_type, self.alphas, self.con_coeffs[:,i])
            for i in xrange(self.con_coeffs.shape[1])
        ]


class GOBasis(object):
    """This class describes basis sets applied to a certain molecular structure."""
    def __init__(self, shell_map, nprims, shell_types, alphas, con_coeffs):
        """
           **Arguments:**

           shell_map
                An array with the center index for each shell.
                shape = (nshell,)

           nprims
                The number of primitives in each shell.
                shape = (nshell,)

           shell_types
                An array with contraction types: 0 = S, 1 = P, 2 = Cartesian D,
                3 = Cartesian F, ..., -2 = pure D, -3 = pure F, ...
                shape = (nshell,)

           alphas
                The exponents of the primitives in one shell.
                shape = (sum(nprims),)

           con_coeffs
                The contraction coefficients of the primitives for each
                contraction in a contiguous array. The coefficients are ordered
                according to the shells. Within each shell, the coefficients are
                grouped per exponent.
                shape = (sum(nprims),)

           Convention for basis functions of a given contraction type:

           The order of the pure shells is based on the order of real spherical.
           The functions are sorted from low to high magnetic quantum number,
           with cosine-like functions before the sine-like functions. The order
           of functions in a Cartesian shell is alhpabetic. Some examples:

           shell_type = 0, S:
             0 -> 1
           shell_type = 1, P:
             0 -> x
             1 -> y
             2 -> z
           shell_type = 2, Cartesian D:
             0 -> xx
             1 -> xy
             2 -> xz
             3 -> yy
             4 -> yz
             5 -> zz
           shell_type = 3, Cartesian F:
             0 -> xxx
             1 -> xxy
             2 -> xxz
             3 -> xyy
             4 -> xyz
             5 -> xzz
             6 -> yyy
             7 -> yyz
             8 -> yzz
             9 -> zzz
           shell_type = -1, not allowed
           shell_type = -2, pure D:
             0 -> zz
             1 -> xz
             2 -> yz
             3 -> xx-yy
             4 -> xy
           shell_type = -3, pure F:
             0 -> zzz
             1 -> xzz
             2 -> yzz
             3 -> xxz-yyz
             4 -> xyz
             5 -> xxx-3xyy
             6 -> 3xxy-yyy
        """
        # All fields are stored as internal parameters. Once they are set,
        # they are no supposed to be modified.
        self._shell_map = shell_map
        self._nprims = nprims
        self._shell_types = shell_types
        self._alphas = alphas
        self._con_coeffs = con_coeffs
        # derived property, read only
        self._nbasis = sum(get_shell_nbasis(shell_type) for shell_type in shell_types)

    def get_nshell(self):
        return len(self._shell_map)

    nshell = property(get_nshell)

    def get_nbasis(self):
        return self._nbasis

    nbasis = property(get_nbasis)

    def compute_overlap(self, centers, overlap):
        from horton.cext import compute_gobasis_overlap
        compute_gobasis_overlap(
            centers, self._shell_map, self._nprims, self._shell_types,
            self._alphas, self._con_coeffs, overlap._array
        )
