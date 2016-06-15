#!/usr/bin/env python
# -*- coding: utf-8 -*-
# HORTON: Helpful Open-source Research TOol for N-fermion systems.
# Copyright (C) 2011-2016 The HORTON Development Team
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
"""Pre-computation and code generation for all sorts of moments"""


from horton import get_cartesian_powers, get_ncart


def write_cart_poly_code(f):
    cps = get_cartesian_powers(7)
    cps = [tuple(row) for row in cps[1:]] # chop of the s function

    print >>f, "    // Shell l=0"
    print >>f, "    if (lmax <= 0) return -1;"
    print >>f


    last_shell = 0
    begin = -1
    end = 0
    for px, py, pz in cps:
        shell = px+py+pz
        if shell > last_shell:
            print >>f, "    // Shell l=%i" % shell
            begin = end
            last_shell = shell

        if shell > 1:
            if shell > 2:
                if px != 0:
                    prev_row = (px-1, py, pz)
                    symbol = 'output[0]'
                elif py != 0:
                    prev_row = (px, py-1, pz)
                    symbol = 'output[1]'
                else:
                    prev_row = (px, py, pz-1)
                    symbol = 'output[2]'
                formula = '%s*output[%i]' % (symbol, cps.index(prev_row))
            else:
                symbols = ["output[0]"]*px + ["output[1]"]*py + ["output[2]"]*pz
                formula = "*".join(symbols)
            print >>f, "    output[%i] = %s;" % (end, formula)
        end += 1

        if end >= begin+get_ncart(shell):
            print >>f, "    if (lmax <= %i) return %i;" % (shell, begin)
            print >>f

    print >> f, '    throw std::domain_error("Encountered lmax > 7.");'


def write_cart_rotation_code(f):
    import numpy as np
    from sympy import Symbol, Wild, Add, Mul, Integer


    R = [Symbol('R_%d' % i) for i in xrange(9)]
    x = Symbol('x')
    y = Symbol('y')
    z = Symbol('z')

    lmax = 4
    cartesian_powers = [tuple(row) for row in get_cartesian_powers(lmax)]
    print cartesian_powers
    ncart = len(cartesian_powers)


    transforms = []

    def get_power(term, base):
        p = Wild("p")
        c = Wild("c")
        d = term.match(c*base**p)
        if d is None:
            return 0
        else:
            return int(d[p])

    icart = 0
    for ishell in xrange(lmax+1):
        shell_rules = []
        transforms.append(shell_rules)
        for ifn in xrange(((ishell+1)*(ishell+2))/2):
            px0, py0, pz0 = cartesian_powers[icart]
            rules = []
            shell_rules.append(rules)
            poly = (
                 (R[0]*x + R[1]*y + R[2]*z)**px0
                *(R[3]*x + R[4]*y + R[5]*z)**py0
                *(R[6]*x + R[7]*y + R[8]*z)**pz0
            ).expand()
            print poly

            if isinstance(poly, Add):
                terms = poly.args
            else:
                terms = [poly]

            for term in terms:
                px1 = get_power(term, x)
                py1 = get_power(term, y)
                pz1 = get_power(term, z)
                prs = [get_power(term, R[i]) for i in xrange(9)]
                col = cartesian_powers.index((px1, py1, pz1)) - cartesian_powers.index((px1+py1+pz1, 0, 0))
                assert col >= 0
                coeff = 1
                if isinstance(term, Mul):
                    for factor in term.args:
                        if isinstance(factor, Integer):
                            coeff *= int(factor)
                rule = (col, coeff, prs)
                rules.append(rule)

            icart += 1


    def format_rule(rule):
        l = [rule[0], rule[1]]
        prs = rule[2]
        for j in xrange(9):
            for p in xrange(prs[j]):
                l.append(j)
            #if prs[j] > 0:
            #    l.append(j)
            #    l.append(prs[j])
        return '[%s]' % (', '.join('%2i' % i for i in l))


    print >> f, 'cartesian_transforms = ['
    for shell_rules in transforms:
        print >> f, '  ['
        for rules in shell_rules:
            if len(rules) == 1:
                print >> f, '    [%s],' % format_rule(rules[0])
            else:
                rules = sorted(rules)
                print >> f, '    [%s,' % format_rule(rules[0])
                for rule in rules[1:-1]:
                    print >> f, '     %s,' % format_rule(rule)
                print >> f, '     %s],' % format_rule(rules[-1])
        print >> f, '  ],'
    print >> f, ']'


if __name__ == "__main__":
    #with file("cart_poly.inc.cpp", "w") as f:
    #    write_cart_poly_code(f)
    with file("cart_rotate.inc.py", "w") as f:
        write_cart_rotation_code(f)
