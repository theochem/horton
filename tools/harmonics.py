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
"""Tools to define the relations between pure and Cartesian gaussian functions.

   This module uses sympy for the symbolic manipulations.
"""


from sympy import sqrt, Symbol, pi, S, Add, cos, sin, Wild
try:
    from sympy.mpmath import fac, fac2, binomial
except ImportError:
    from mpmath import fac, fac2, binomial

import numpy as np


def get_cart_dof(l):
    return ((l+1)*(l+2))/2


def get_pure_dof(l):
    return 2*l+1


def get_cart_orb_norm(alpha, nx, ny, nz):
    """Returns the norm of a Cartesian gaussian funcion"""
    return sqrt(
        int(fac2(2*nx-1)*fac2(2*ny-1)*fac2(2*nz-1))
        /(2*alpha/pi)**(3/S(2))
        /(4*alpha)**(nx+ny+nz)
    )


def get_pure_orb_norm(alpha, l):
    """Returns the norm of a pure gaussian funcion"""
    return sqrt(
        int(fac2(2*l-1))
        /(2*alpha/pi)**(3/S(2))
        /(4*alpha)**l
    )


def iter_cartesian_powers(order):
    """Iterate over Cartesian powers in `alphabetical` order"""
    for nx in xrange(order,-1,-1):
        for ny in xrange(order-nx,-1,-1):
            yield nx, ny, order-nx-ny


def get_cart_orb_polys(shell, alpha, xyz):
    """Returns a list of normalized Cartesian polynomials used in Gaussian basis functions

       **Arguments:**

       shell
            positive shell parameter (0=s, 1=p, 2=d, ...)

       alpha
            The exponent of the Gaussian. This is used for the correct
            normalization.

       xyz
            A tuple with Sympy x, y, z symbols.
    """
    x, y, z = xyz
    result = []
    for nx, ny, nz in iter_cartesian_powers(shell):
        result.append(x**nx*y**n*z**nz/get_cart_orb_norm(alpha, nx, ny, nz))
    return result


def get_pure_orb_polys(shell, alpha, xyz):
    """Returns a list of normalized pure polynomials used in Gaussian basis functions

       **Arguments:**

       shell
            positive shell parameter (0=s, 1=p, 2=d, ...)

       alpha
            The exponent of the Gaussian. This is used for the correct
            normalization.

       xyz
            A tuple with Sympy x, y, z symbols.
    """
    norm = get_pure_orb_norm(alpha, shell)
    x, y, z = xyz
    polys = get_solid_harmonics(shell, sqrt(x*x+y*y+z*z), x, y, z)
    result = []
    for poly in polys:
        result.append((poly/norm).expand())
    return result


def get_solid_harmonics(l, r, x, y, z):
    def get_pi(l, m):
        terms = []
        for k in xrange((l-m)/2+1):
            factor = (-1)**k*int(binomial(l,k)*binomial(2*l-2*k,l)*fac(l-2*k))/S(2**l*int(fac(l-2*k-m)))
            terms.append(factor*r**(2*k)*z**(l-2*k-m))
        return Add(*terms)

    def get_a(m):
        terms = []
        for p in xrange(m+1):
            factor = int(binomial(m, p))*cos((m-p)*pi/S(2))
            terms.append(factor*x**p*y**(m-p))
        return Add(*terms)

    def get_b(m):
        terms = []
        for p in xrange(m+1):
            factor = int(binomial(m, p))*sin((m-p)*pi/S(2))
            terms.append(factor*x**p*y**(m-p))
        return Add(*terms)

    result = []
    for m in xrange(l+1):
        if m == 0:
            result.append(get_pi(l, m)*get_a(m))
        else:
            factor = sqrt(2*int(fac(l-m))/S(int(fac(l+m))))
            result.append(factor*get_pi(l, m)*get_a(m))
            result.append(factor*get_pi(l, m)*get_b(m))

    return result


def get_poly_conversion(shell):
    # conversion from normalised Cartesian to normalized pure
    alpha = Symbol("alpha")
    x = Symbol("x")
    y = Symbol("y")
    z = Symbol("z")
    xyz = (x,y,z)

    def get_power(term, base):
        p = Wild("p")
        c = Wild("c")
        d = term.match(c*base**p)
        if d is None:
            return 0
        else:
            return int(d[p])

    num_dof_in = get_cart_dof(shell)
    num_dof_out = get_pure_dof(shell)
    lcs = np.zeros((num_dof_out, num_dof_in), dtype=object)

    for i_out, poly in enumerate(get_pure_orb_polys(shell, alpha, xyz)):
        poly = poly.expand()
        if isinstance(poly, Add):
            terms = poly.args
        else:
            terms = [poly]
        coeffs = {}
        for term in terms:
            px = get_power(term, x)
            py = get_power(term, y)
            pz = get_power(term, z)
            key = (px, py, pz)
            value = term.subs({x:1,y:1,z:1})
            coeffs[key] = value
        for i_in, key in enumerate(iter_cartesian_powers(shell)):
            cart_wfn_norm = get_cart_orb_norm(alpha, key[0], key[1], key[2])
            lcs[i_out,i_in] = coeffs.get(key, 0)*cart_wfn_norm
    return lcs


def iter_labels(l):
    m = 0
    yield 'C_%i^%i' % (l,m)
    for counter in xrange(l):
        m += 1
        yield 'C_%i^%i' % (l,m)
        yield 'S_%i^%i' % (l,m)


def run_print_solid_harmonics():
    r = Symbol("r")
    x = Symbol("x")
    y = Symbol("y")
    z = Symbol("z")
    xyz = (x, y, z)

    print
    print '.. math::'
    for shell in xrange(5):
        il = iter_labels(shell)
        for poly in get_solid_harmonics(shell, r, x, y, z):
            label = il.next()
            print '    %s(x,y,z) & = %s \\\\' % (label, latex(poly))#.replace('$',''))


def run_print_transformations_latex():
    x = Symbol("x")
    y = Symbol("y")
    z = Symbol("z")
    xyz = (x, y, z)

    print
    print '.. math::'
    for shell in xrange(5):
        lcs = get_poly_conversion(shell)
        npure, ncart = lcs.shape
        print '   ', r'\left(\begin{array}{c}'
        print '   ', r' \\ '.join(['X(%s)' % label for label in iter_labels(shell)])
        print '   ', r'\end{array}\right)'
        print '    &='
        print '   ', r'\left(\begin{array}{%s}' % ('c'*ncart)
        for ipure in xrange(npure):
            print '   ', r' & '.join([latex(lcs[ipure, icart]) for icart in xrange(ncart)]), r'\\'
        print '   ', r'\end{array}\right)'
        print '   ', r'\left(\begin{array}{c}'
        els = []
        for nx, ny, nz in iter_cartesian_powers(shell):
            spoly = 'x'*nx+'y'*ny+'z'*nz
            if spoly == '': spoly = '1'
            els.append('X(%s)' % spoly)
        print '   ', r' \\ '.join(els)
        print '   ', r'\end{array}\right)'
        print r'    \\'


def strip_zero(s):
    while len(s) > 3 and s[-1] == '0' and s[-2] != '.':
        s = s[:-1]
    return s


def run_print_transformations_c():
    x = Symbol("x")
    y = Symbol("y")
    z = Symbol("z")
    xyz = (x, y, z)

    nshell = 10
    sizes = []

    print
    for shell in xrange(nshell):
        print 'const static type_sparse_el cptf%i[] = {' % shell
        lcs = get_poly_conversion(shell)
        npure, ncart = lcs.shape
        size = 0
        for ipure in xrange(npure):
            for icart in xrange(ncart):
                val = lcs[ipure,icart].evalf(20)
                if val != 0:
                    val = strip_zero(str(val))
                    print '    {%2i, %2i, %s},' % (ipure, icart, val)
                    size +=1
        print '};'
        sizes.append(size)
    print 'const static type_sparse_tf cptf[MAX_CON_TYPE+1] = {'
    print '    %s' % (', '.join('{%i, cptf%i}' % (sizes[shell], shell) for shell in xrange(nshell)))
    print '};'


def run_print_transformations_python():
    x = Symbol("x")
    y = Symbol("y")
    z = Symbol("z")
    xyz = (x, y, z)

    nshell = 5

    print
    for shell in xrange(nshell):
        print 'tf = np.array(['
        lcs = get_poly_conversion(shell)
        npure, ncart = lcs.shape
        size = 0
        for ipure in xrange(npure):
            vals = []
            for icart in xrange(ncart):
                val = lcs[ipure,icart].evalf(20)
                val = strip_zero(str(val))
                vals.append(val)
            print '    [%s],' % (', '.join(vals))
        print '])'


if __name__ == "__main__":
    #run_print_solid_harmonics()
    #run_print_transformations_latex()
    #run_print_transformations_c()
    run_print_transformations_python()
