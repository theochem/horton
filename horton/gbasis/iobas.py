# -*- coding: utf-8 -*-
# HORTON: Helpful Open-source Research TOol for N-fermion systems.
# Copyright (C) 2011-2015 The HORTON Development Team
#
# This file is part of HORTON.
#
# HORTON is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
#
# HORTON is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, see <http://www.gnu.org/licenses/>
#
#--
'''Input/Output routines for gaussian basis sets'''

import numpy as np

from horton.periodic import periodic


__all__ = [
    'str_to_shell_types', 'shell_type_to_str', 'fortran_float',
    'load_basis_atom_map_nwchem', 'load_basis_atom_map_gbs',
    'write_gbs'
]


def str_to_shell_types(s, pure=False):
    """Convert a string into a list of contraction types"""
    if pure:
        d = {'s': 0, 'p': 1, 'd': -2, 'f': -3, 'g': -4, 'h': -5, 'i': -6}
    else:
        d = {'s': 0, 'p': 1, 'd': 2, 'f': 3, 'g': 4, 'h': 5, 'i': 6}
    return [d[c] for c in s.lower()]


def shell_type_to_str(shell_type):
    """Convert a shell type into a character"""
    return {0: 's', 1: 'p', 2: 'd', 3: 'f', 4: 'g', 5: 'h', 6: 'i'}[abs(shell_type)]


def fortran_float(s):
    '''Convert a string to a float. Works also with D before the mantissa'''
    return float(s.replace('D', 'E').replace('d', 'e'))


def load_basis_atom_map_nwchem(filename):
    '''Load the basis set family from an NWChem file.'''
    from horton.gbasis.gobasis import GOBasisAtom, GOBasisContraction

    f = open(filename)
    basis_atom_map = {}
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
            ba = basis_atom_map.get(n)
            shell_types = str_to_shell_types(words[1])
            bcs = [GOBasisContraction(shell_type, [], []) for shell_type in shell_types]
            if ba is None:
                ba = GOBasisAtom(bcs)
                basis_atom_map[n] = ba
            else:
                ba.bcs.extend(bcs)
        else:
            # An extra primitive for the current contraction(s).
            exponent = fortran_float(words[0])
            coeffs = [fortran_float(w) for w in words[1:]]
            for i, bc in enumerate(bcs):
                bc.alphas.append(exponent)
                bc.con_coeffs.append(coeffs[i::len(bcs)])
    f.close()
    return basis_atom_map

def load_basis_atom_map_gbs(filename):
    """ Loads the basis set family from a GBS file

    """
    from horton.gbasis.gobasis import GOBasisAtom, GOBasisContraction

    basis_atom_map = {}
    bc = None # The current contraction being loaded
    cur_atom = None
    cur_shell_types = None
    with open(filename, 'r') as f:
        for line in f:
            # strip of comments and white space
            line = line[:line.find('!')].strip()
            if len(line) == 0:
                continue
            if line == '****':
                continue
            words = line.split()
            # if first word is the atomic symbol
            if words[0].isalpha() and len(words) == 2:
                cur_atom = words[0]
            # if first word is the angular momentum
            elif words[0].isalpha() and len(words) == 3:
                # A new contraction begins, maybe even a new atom.
                n = periodic[cur_atom].number
                cur_shell_types = str_to_shell_types(words[0])
                empty_contr = [GOBasisContraction(shell_type, [], []) for shell_type in cur_shell_types]
                try:
                    basis_atom_map[n].bcs.extend(empty_contr)
                except KeyError:
                    basis_atom_map[n] = GOBasisAtom(empty_contr)
            else:
                # An extra primitive for the current contraction(s).
                exponent = fortran_float(words[0])
                coeffs = [fortran_float(w) for w in words[1:]]
                for i, bc in enumerate(empty_contr):
                    bc.alphas.append(exponent)
                    bc.con_coeffs.append(coeffs[i::len(cur_shell_types)])
    return basis_atom_map

def write_gbs(filename, gaussian_basis):
    """ Writes gaussian basis file from the basis object in HORTON

    Parameters
    ----------
    filename: str
        File name of the new gbs file
    gaussian_basis: GOBasisFamily
        Basis set object that contains the name and the contraction information of
        the basis set

    Raises
    ------
    AssertionError
        If gaussian_basis is not loaded

    """
    with open(filename, 'w') as f:
        f.write('!Basis set, {0}, generated using HORTON\n\n'.format(gaussian_basis.name))
        f.write('****\n')
        basis_atom_map = gaussian_basis.basis_atom_map
        assert not basis_atom_map is None, ('GOBasisFamily object was not loaded.'
                                            ' Fix with `gaussian_basis.load()` before'
                                            ' passing the gaussian basis object')
        for atom in sorted(basis_atom_map.keys()):
            f.write('{0:<6}0\n'.format(periodic[atom].symbol))
            contractions = basis_atom_map[atom].bcs
            for contraction in contractions:
                exponents = contraction.alphas.reshape(contraction.alphas.size, 1)
                con_coeffs = contraction.con_coeffs
                if len(contraction.con_coeffs.shape) == 1:
                    con_coeffs = contraction.con_coeffs.reshape(contraction.con_coeffs.size, 1)
                con_numbers = np.hstack((exponents, con_coeffs))
                f.write('{0:<4}{1:<4}1.00\n'.format(shell_type_to_str(contraction.shell_type).upper(),
                                                    exponents.size))
                for con_number in con_numbers:
                    f.write(('{:>17}'*con_number.size).format(*con_number))
                    f.write('\n')
            f.write('****\n')
