# -*- coding: utf-8 -*-
# Horton is a Density Functional Theory program.
# Copyright (C) 2011-2013 Toon Verstraelen <Toon.Verstraelen@UGent.be>
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


from horton.periodic import periodic


__all__ = ['str_to_shell_types', 'load_basis_atom_map_nwchem']


def str_to_shell_types(s, pure=False):
    """Convert a string into a list of contraction types"""
    if pure:
        d = {'s': 0, 'p': 1, 'd': -2, 'f': -3, 'g': -4, 'h': -5, 'i': -6}
    else:
        d = {'s': 0, 'p': 1, 'd': 2, 'f': 3, 'g': 4, 'h': 5, 'i': 6}
    return [d[c] for c in s.lower()]


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
            exponent = float(words[0])
            coeffs = [float(w) for w in words[1:]]
            for i, bc in enumerate(bcs):
                bc.alphas.append(exponent)
                bc.con_coeffs.append(coeffs[i::len(bcs)])
    f.close()
    return basis_atom_map
