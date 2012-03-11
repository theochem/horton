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


__all__ = [
    'GOBasisDesc', 'str_to_con_types', 'GOBasisFamily', 'go_basis_families',
    'GOBasisAtom', 'GOBasisShell', 'get_con_nbasis', 'GOBasis'
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
        num_exponents = []
        num_contractions = []
        con_types = []
        con_coeffs = []
        exponents = []

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
                    raise ValueError('The basis family %s does not contain element %n.' % (basis_x, n))
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
            basis_atom.extend(i, shell_map, num_exponents, num_contractions, con_types, exponents, con_coeffs, self.pure)

        # Turn the lists into arrays
        shell_map = np.array(shell_map)
        num_exponents = np.array(num_exponents)
        num_contractions = np.array(num_contractions)
        con_types = np.array(con_types)
        con_coeffs = np.array(con_coeffs)
        exponents = np.array(exponents)

        # Return the Gaussian basis object.
        return GOBasis(shell_map, num_exponents, num_contractions, con_types, exponents, con_coeffs)


def str_to_con_types(s):
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
        f = open(self.filename)
        self.basis_atom_map = {}
        bc = None # The current basis contraction being loaded
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
                con_types = str_to_con_types(words[1])
                bc = GOBasisShell(con_types, [], [])
                if ba is None:
                    ba = GOBasisAtom([bc])
                    self.basis_atom_map[n] = ba
                else:
                    ba.contractions.append(bc)
            else:
                # An extra primitive for the current contraction.
                exponent = float(words[0])
                coeffs = [float(w) for w in words[1:]]
                bc.exponents.append(exponent)
                bc.con_coeffs.extend(coeffs)
        f.close()


go_basis_families = [
    GOBasisFamily('STO-3G', filename=context.get_fn('basis/sto-3g.nwchem')),
    GOBasisFamily('3-21G', filename=context.get_fn('basis/3-21g.nwchem')),
]
go_basis_families = dict((bf.name.lower(), bf) for bf in go_basis_families)


class GOBasisAtom(object):
    """Description of an atomic basis set with generalized contractions."""
    def __init__(self, contractions):
        self.contractions = contractions

    def extend(self, i, shell_map, num_exponents, num_contractions, con_types, exponents, con_coeffs, pure=True):
        """Extend the lists with basis functions for this atom.

           **Arguments:**

           i
                The index of the center of this atom.

           shell_map, num_exponents, num_contractions, con_types, con_coeffs, exponents
                These are all output arguments that must be extended with the
                basis set parameters for this atom. The meaning of each argument
                is defined in the documentation for constructor of the GOBasis
                class.

           **Optional argument:**

           pure
                By default pure basis functions are used. Set this to false to
                switch to Cartesian basis functions.
        """
        for bc in self.contractions:
            shell_map.append(i)
            num_exponents.append(len(bc.exponents))
            num_contractions.append(len(bc.con_types))
            if pure:
                con_types.extend(bc.con_types)
            else:
                for con_type in bc.con_types:
                    # swap sign for d and higher.
                    if con_type > 1:
                        con_type *= -1
                    con_types.append(con_type)
            exponents.extend(bc.exponents)
            for coeffs in bc.con_coeffs:
                con_coeffs.append(coeffs)


class GOBasisShell(object):
    """Description of a generalized contraction for a given shell"""
    def __init__(self, con_types, exponents, con_coeffs):
        """
           **Arguments:**

           con_types
                a list of integers, the angular quantum numbers of each
                contraction in the shell. (0=s, 1=p, 2=d, ...)

           exponents
                The exponents of the primitives in each contraction. These
                are shared between the different contractions in the
                shell. (This is the idea of the generalized contraction.)

           con_coeffs
                A 2D array (or list of rows) with all the linear coefficients
                of the contractions. Each row corresponds to one exponent.
                Each column corresponds to one contraction.
        """
        self.con_types = con_types
        self.exponents = exponents
        self.con_coeffs = con_coeffs


def get_con_nbasis(con_type):
    """Return the number basis function for a given contraction type."""
    if con_type > 0:
        return (con_type+1)*(con_type+2)/2
    elif con_type == -1:
        raise ValueError
    else:
        return -2*con_type+1


class GOBasis(object):
    """This class describes basis sets applied to a certain molecular structure."""
    def __init__(self, shell_map, num_exponents, num_contractions, con_types, exponents, con_coeffs):
        """
           **Arguments:**

           shell_map
                An array with the center index for each shell.

           num_exponents
                The number of exponents in each shell.

           num_contractions
                The number of contractions in each shell. This is used to
                implement optimized general contractions.

           con_types
                An array with contraction types: 0 = S, 1 = P, 2 = Cartesian D,
                3 = Cartesian F, ..., -2 = pure D, -3 = pure F, ...
                One contraction type is present for each contraction in each
                shell. The so-called SP type is implemented as a shell
                with two contractions, one of type S and one of type P.

           exponents
                The exponents of the primitives in one shell.

           con_coeffs
                The contraction coefficients of the primitives for each
                contraction in a contiguous array. The coefficients are ordered
                according to the shells. Within each shell, the coefficients are
                grouped per exponent.

           The number of primitives in shell i is num_exponents[i]*num_contractions[i].

           Convention for basis functions of a given contraction type:

           The order of the pure shells is based on the order of real spherical
           harmonics: http://en.wikipedia.org/wiki/Table_of_spherical_harmonics
           First the +- linear combination of highest angular momentum, then
           the ++ combination of highest angular momentum, keep repeating and
           finally take angular momention zero (without making a linear
           combination). The order of the Cartesian shells is sorted
           alhpabetically. The SP shell type is S first, then P. Some examples:

           con_type=0, S:
             0 -> 1
           con_type=1, P:
             0 -> x
             1 -> y
             2 -> z
           con_type=2, Cartesian D:
             0 -> xx
             1 -> xy
             2 -> xz
             3 -> yy
             4 -> yz
             5 -> zz
           con_type=3, Cartesian F:
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
           con_type=-1, SP:
             0 -> 1
             1 -> x
             2 -> y
             3 -> z
           con_type=-2, pure D:
             0 -> zz
             1 -> yz
             2 -> xz
             3 -> xx-yy
             4 -> xy
           con_type=-3, pure F:
             6 -> zzz
             5 -> yzz
             4 -> xzz
             3 -> xxz-yyz
             2 -> xyz
             1 -> 3xxy-yyy
             0 -> xxx-3xyy
        """
        # All fields are stored as internal parameters. Once they are set,
        # they are no supposed to be modified.
        self._shell_map = shell_map
        self._num_exponents = num_exponents
        self._num_contractions = num_contractions
        self._con_types = con_types
        self._exponents = exponents
        self._con_coeffs = con_coeffs
        # derived property, read only
        self._nbasis = sum(get_con_nbasis(con_type) for con_type in con_types)

    def get_nshell(self):
        return len(self._shell_map)

    nshell = property(get_nshell)

    def get_nbasis(self):
        return self._nbasis

    nbasis = property(get_nbasis)

    @classmethod
    def from_fchk(cls, filename):
        from horton.io import FCHKFile
        fchk = FCHKFile(filename, [
            "Shell types", "Shell to atom map", "Shell to atom map",
            "Number of primitives per shell", "Primitive exponents",
            "Contraction coefficients", "P(S=P) Contraction coefficients",
        ])
        shell_types = fchk.fields["Shell types"]
        shell_map = fchk.fields["Shell to atom map"] - 1
        num_exponents = fchk.fields["Number of primitives per shell"]
        exponents = fchk.fields["Primitive exponents"]
        ccoeffs_level1 = fchk.fields["Contraction coefficients"]
        ccoeffs_level2 = fchk.fields.get("P(S=P) Contraction coefficients")

        num_contractions = []
        con_types = []
        con_coeffs = []
        counter = 0
        for i, n in enumerate(num_exponents):
            if shell_types[i] == -1:
                # Special treatment for SP shell type
                num_contractions.append(2)
                con_types.append(0)
                con_types.append(1)
                tmp = np.array([
                    ccoeffs_level1[counter:counter+n],
                    ccoeffs_level2[counter:counter+n]
                ])
                con_coeffs.append(tmp.transpose().ravel())
            else:
                num_contractions.append(1)
                con_types.append(shell_types[i])
                con_coeffs.append(ccoeffs_level1[counter:counter+n])
            counter += n
        num_contractions = np.array(num_contractions)
        con_types = np.array(con_types)
        con_coeffs = np.concatenate(con_coeffs)


        result = cls(shell_map, num_exponents, num_contractions, con_types, exponents, con_coeffs)

        # permutation of the basis functions
        g03_reordering = {
          -6: np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]),
          -5: np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]),
          -4: np.array([0, 1, 2, 3, 4, 5, 6, 7, 8]),
          -3: np.array([0, 1, 2, 3, 4, 5, 6]),
          -2: np.array([0, 1, 2, 3, 4]),
          -1: np.array([0, 1, 2, 3]),
           0: np.array([0]),
           1: np.array([0, 1, 2]),
           2: np.array([0, 3, 4, 1, 5, 2]),
           3: np.array([0, 4, 5, 3, 9, 6, 1, 8, 7, 2]),
           4: np.arange(15)[::-1],
           5: np.arange(21)[::-1],
           6: np.arange(28)[::-1],
        }
        offset = 0
        permutation = []
        for shell_type in shell_types:
            permutation.extend(g03_reordering[shell_type]+len(permutation))
        result.g03_permutation = np.array(permutation, dtype=int)

        return result

    def compute_overlap(self, centers, overlap):
        from horton.ext import compute_gobasis_overlap
        compute_gobasis_overlap(
            centers, self._shell_map, self._num_exponents,
            self._num_contractions, self._con_types, self._exponents,
            self._con_coeffs, overlap._array
        )
