# -*- coding: utf-8 -*-
# HORTON: Helpful Open-source Research TOol for N-fermion systems.
# Copyright (C) 2011-2017 The HORTON Development Team
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
# --
"""Input/Output routines for gaussian basis sets"""

import numpy as np

from .periodic import num2sym, sym2num

__all__ = [
    'str_to_shell_types', 'shell_type_to_str', 'fortran_float',
    'load_basis_atom_map_nwchem', 'load_basis_atom_map_gbs',
    'dump_basis_atom_map_gbs'
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
    """Convert a string to a float. Works also with D before the mantissa"""
    return float(s.replace('D', 'E').replace('d', 'e'))


def load_basis_atom_map_nwchem(filename):
    """Load the basis set family from an NWChem file."""
    from .gobasis import GOBasisAtom, GOBasisContraction

    f = open(filename)
    basis_atom_map = {}
    for line in f:
        # Strip off comments and white space.
        line = line[:line.find('#')].strip()
        if len(line) == 0 or line.startswith('BASIS'):
            continue
        if line == 'END':
            break
        words = line.split()
        if words[0].isalpha():
            # A new contraction begins, maybe even a new atom.
            n = sym2num[words[0]]
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
    """Load the basis set family from a GBS file."""
    from .gobasis import GOBasisAtom, GOBasisContraction

    basis_atom_map = {}
    cur_atom = None
    cur_shell_types = None
    with open(filename, 'r') as f:
        for line in f:
            # Strip off comments and white space.
            line = line[:line.find('!')].strip()
            if len(line) == 0 or line == '****':
                continue
            words = line.split()
            # if first word is the atomic symbol
            if words[0].isalpha() and len(words) == 2:
                cur_atom = words[0]
            # if first word is the angular momentum
            elif words[0].isalpha() and len(words) == 3:
                # A new contraction begins, maybe even a new atom.
                n = sym2num[cur_atom]
                cur_shell_types = str_to_shell_types(words[0])
                empty_contr = [GOBasisContraction(shell_type, [], [])
                               for shell_type in cur_shell_types]
                # Try to get the atom and add empty contraction, or create new atom.
                goba = basis_atom_map.get(n)
                if goba is None:
                    basis_atom_map[n] = GOBasisAtom(empty_contr)
                else:
                    goba.bcs.extend(empty_contr)
            else:
                # An extra primitive for the current contraction(s).
                exponent = fortran_float(words[0])
                coeffs = [fortran_float(w) for w in words[1:]]
                for i, bc in enumerate(empty_contr):  # FIXME: empty_contr ref before assignment
                    bc.alphas.append(exponent)
                    bc.con_coeffs.append(coeffs[i::len(cur_shell_types)])
    return basis_atom_map


def dump_basis_atom_map_gbs(filename, name, basis_atom_map):
    """Write gaussian basis file from the basis object in HORTON.

    Parameters
    ----------
    filename: str
        File name of the new gbs file
    name : str
        Name of the basis set to mention in the comments of the written file.
    basis_atom_map: dict
        Keys are atomic numbers, values are GOBasisAtom objects.
    """
    with open(filename, 'w') as f:
        f.write('!Basis set, {0}, generated using HORTON\n\n'.format(name))
        f.write('****\n')
        for atom, gobatom in sorted(basis_atom_map.items()):
            f.write('{0:<6}0\n'.format(num2sym[atom]))
            contractions = gobatom.bcs
            for contraction in contractions:
                exponents = contraction.alphas.reshape(-1, 1)
                con_coeffs = contraction.con_coeffs
                if con_coeffs.ndim == 1:
                    con_coeffs = contraction.con_coeffs.reshape(-1, 1)
                con_numbers = np.hstack((exponents, con_coeffs))
                f.write('{0:<4}{1:<4}1.00\n'.format(
                    shell_type_to_str(contraction.shell_type).upper(), exponents.size))
                for con_row in con_numbers:
                    f.write(('{:>17} ' * con_row.size).format(*con_row))
                    f.write('\n')
            f.write('****\n')
