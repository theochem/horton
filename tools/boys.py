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
"""Pre-computes function values of the Boys function to full precision.
"""


from sympy import *
try:
    from sympy.mpmath.libmp import to_str
except ImportError:
    from mpmath.libmp import to_str



def boys_mpmath(m, t):
    if t == 0:
        return 1/(2*m+1)
    else:
        eps = m+Rational(1,2)
        return lowergamma(eps, t)/(2*(t**eps))


def boys_function_tail(m, t):
    result = 1;
    for i in xrange(1, m+1):
        result *= (2*i-1)/(2*t)
    return result*sqrt(pi/t)/2


def print_boys_inc(maxm, resolution):
    # Computation of all the boys function data
    boys_fn_data = []
    for m in xrange(maxm+1):
        l = []
        boys_fn_data.append(l)
        m = sympify(m)
        counter = 0
        while True:
            t = Rational(counter, resolution)
            f = boys_mpmath(m, t)
            f_numer = f.evalf(30, maxn=1000, strict=True)
            f_str = to_str(f_numer._mpf_, 16, strip_zeros=False, min_fixed=0, max_fixed=0)
            l.append(f_str)
            ft = boys_function_tail(m, t)
            ft_numer = ft.evalf(30, maxn=1000, strict=True)
            counter += 1
            error = abs(float(ft) - float(f))
            if counter % 100 == 0 or error == 0:
                print '%2i %5i % 20.12e % 20.12e' % (m, counter, f, error)
            if error == 0:
                break

    f = file('boys_inc.cpp', 'w')
    print >> f
    print >> f
    print >> f, '#define BOYS_RESOLUTION', resolution
    print >> f, '#define BOYS_MAX_DATA', maxm
    print >> f
    for m in xrange(maxm+1):
        print >> f, 'static double boys_fn_data_%i[%i] = {' % (m, len(boys_fn_data[m]))
        print >> f, wrap_strings(boys_fn_data[m], 22, 10)
        print >> f, '};'
        print >> f

    print >> f, 'static double *boys_fn_data[BOYS_MAX_DATA+1] = {'
    print >> f, wrap_strings(['boys_fn_data_%i' % m for m in xrange(maxm+1)], 16, 4)
    print >> f, '};'
    print >> f

    print >> f, 'static long boys_sizes[BOYS_MAX_DATA+1] = {'
    print >> f, wrap_strings([str(len(boys_fn_data[m])) for m in xrange(maxm+1)], 6, 8)
    print >> f, '};'
    print >> f

    f.close()


def wrap_strings(ss, width, ncols):
    result = []
    counter = 0
    for s in ss:
        s = s.rjust(width)
        if counter % ncols == 0:
            result.append('\n    %s,' % s)
        else:
            result.append('%s,' % s)
        counter += 1
    return (''.join(result))[1:-1]


def print_boys_test(maxm):
    l1 = []
    l2 = []

    ms = [S(m) for m in xrange(maxm+1)]
    ts = [S(0), S(10)**(-20), S(10)**(-10), S(10)**(-7), S(10)**(-5),
          S(10)**(-2), S(10)**(-1), 5*S(10)**(-1), S(1), Rational(3,2), S(2),
          S(10), S(10)**2, S(10)**5, S(10)**10, S(10)**20]

    for m in ms:
        for t in ts:
            r = boys_mpmath(m, t).evalf(30, maxn=1000000, strict=True)
            s = to_str(r._mpf_, 20, strip_zeros=False, min_fixed=0, max_fixed=0)
            #print m, t, s
            l1.append('boys_function(%i,%.1e)' % (m, t))
            l2.append(s)


    print 'result = np.array([%s])' % ', '.join(l1)
    print 'check = np.array([%s])' % ', '.join(l2)

if __name__ == '__main__':
    maxm = 4*7
    #print_boys_inc(maxm=maxm+6, resolution=20)
    print_boys_test(maxm=maxm)
