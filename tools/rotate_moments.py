#!/usr/bin/env python
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
"""Pre-computes linear transformation matrices for multipole moments"""


import numpy as np
from sympy import Symbol, Wild, Add, Mul, Integer

from horton import get_cartesian_powers, get_ncart


R = [Symbol('R_%d' % i) for i in xrange(9)]
x = Symbol('x')
y = Symbol('y')
z = Symbol('z')

lmax = 4
cartesian_powers = [tuple(row) for row in get_cartesian_powers(lmax)]
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


for px0, py0, pz0 in cartesian_powers:
    rules = []
    transforms.append(rules)
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
        col = cartesian_powers.index((px1, py1, pz1))
        coeff = 1
        if isinstance(term, Mul):
            for factor in term.args:
                if isinstance(factor, Integer):
                    coeff *= int(factor)
        rule = (col, coeff, prs)
        print rule
        rules.append(rule)

    print


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


print 'cartesian_transforms = ['
for rules in transforms:
    if len(rules) == 1:
        print '    [%s],' % format_rule(rules[0])
    else:
        rules = sorted(rules)
        print '    [%s,' % format_rule(rules[0])
        for rule in rules[1:-1]:
            print '     %s,' % format_rule(rule)
        print '     %s],' % format_rule(rules[-1])
print ']'
