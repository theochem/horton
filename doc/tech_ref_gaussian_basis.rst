..
    : HORTON: Helpful Open-source Research TOol for N-fermion systems.
    : Copyright (C) 2011-2019 The HORTON Development Team
    :
    : This file is part of HORTON.
    :
    : HORTON is free software; you can redistribute it and/or
    : modify it under the terms of the GNU General Public License
    : as published by the Free Software Foundation; either version 3
    : of the License, or (at your option) any later version.
    :
    : HORTON is distributed in the hope that it will be useful,
    : but WITHOUT ANY WARRANTY; without even the implied warranty of
    : MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    : GNU General Public License for more details.
    :
    : You should have received a copy of the GNU General Public License
    : along with this program; if not, see <http://www.gnu.org/licenses/>
    :
    : --

Gaussian basis sets
###################

.. _ref_gaussian_basis_standard_sets:

Standard basis sets
===================

A list of currently supported built-in basis sets is tabulated below.

.. cssclass:: table-striped

.. include:: basis.rst.inc

Note that the basis set names are case-insensitive in HORTON. These basis sets
were taken from the EMSL library (https://bse.pnl.gov/bse/portal). When
publishing results obtained with these basis sets, please cite the following
references [feller1996]_ and [Didier2007]_.


Collected notes on Gaussian basis sets
======================================

Introduction
------------


HORTON supports contracted Gaussian basis functions, which have in general the
following form:

.. math:: b(\mathbf{r}; D_1, \ldots, D_k, P, \alpha_1, \ldots, \alpha_K, \mathbf{r}_A) =
          \sum_{k=1}^K D_k N(\alpha_k, P)
          P(\mathbf{r} - \mathbf{r}_A)
          \exp(-\alpha_k \Vert \mathbf{r} - \mathbf{r}_A \Vert^2)

where :math:`K` is the contraction length, :math:`D_k` is a contraction
coefficient, :math:`N` is a normalization constant, :math:`P` is a Cartesian
polynomial, :math:`\alpha_k` is an exponent and :math:`\mathbf{r}_A` is the
center of the basis function. The summation over :math:`k` is
conventionally called a contraction of `primitive Gaussian basis functions`.
The normalization of each primitive depends on both the polynomial and the
exponent and is defined by the following relation:

.. math:: \int \Bigl\vert N(\alpha_k, P) P(\mathbf{r} - \mathbf{r}_A)
               \exp(-\alpha_k \Vert \mathbf{r} - \mathbf{r}_A \Vert^2)
               \Bigr\vert^2 d\mathbf{r} = 1

Likewise, the contraction coefficients are defined such that the total
contraction satisfies the same normalization condition:

.. math:: \int \Bigl\vert
               b(\mathbf{r}; D_1, \ldots, D_k, P, \alpha_1, \ldots, \alpha_K, \mathbf{r}_A)
               \Bigr\vert^2 d\mathbf{r} = 1

There are two common forms of the polynomial: Cartesian and pure (harmonic) basis
functions. Both types will be defined below, together with some conventions
that are needed for the implementation in HORTON.


Cartesian basis functions
-------------------------


When the polynomial consists of a single term as follows:

.. math:: P(x,y,z) = x^{n_x} y^{n_y} z^{n_z}

with :math:`n_x`, :math:`n_y`, :math:`n_z`, zero or positive integer powers, one
speaks of `Cartesian Gaussian basis functions`. One refers to the sum of the
powers as the angular momentum of the Cartesian Gaussian basis. The
normalization constant is:

.. math:: N(\alpha_k, n_x, n_y, n_z) = \sqrt{\frac
        {(2\alpha_k/\pi)^{3/2} (4\alpha_k)^{n_x+n_y+n_z}}
        {(2n_x-1)!! (2n_y-1)!! (2n_z-1)!!}
        }

In practice one combines all basis functions of a given angular momentum (or
algebraic order). A basis specification typically only mentions the total
angular momentum, and it is assumed that all polynomials of that
order are included in the basis set. The number of basis functions, i.e. the
number of polynomials, for a given angular momentum, :math:`n=n_x+n_y+n_z`, is
:math:`(n+1)(n+2)/2`. For the implementation, one must fix a certain ordering of
these polynomials. In HORTON, the ordering is simply alphabetical.

The first five angular momenta and the corresponding polynomials are listed in
the table below.

====== ================ == ===========
Symbol Angular momentum #  Polynomials
====== ================ == ===========
S      0                1  :math:`1`
P      1                3  :math:`x`, :math:`y`, :math:`z`
D      2                6  :math:`xx`, :math:`xy`, :math:`xz`, :math:`yy`, :math:`yz`, :math:`zz`
F      3                10 :math:`xxx`, :math:`xxy`, :math:`xxz`, :math:`xyy`, :math:`xyz`, :math:`xzz`, :math:`yyy`, :math:`yyz`, :math:`yzz`, :math:`zzz`
G      4                15 :math:`xxxx`, :math:`xxxy`, :math:`xxxz`, :math:`xxyy`, :math:`xxyz`, :math:`xxzz`, :math:`xyyy`, :math:`xyyz`, :math:`xyzz`, :math:`xzzz`, :math:`yyyy`, :math:`yyyz`, :math:`yyzz`, :math:`yzzz`, :math:`zzzz`
====== ================ == ===========


Pure or harmonic basis functions
--------------------------------

When the polynomial is a real regular solid harmonic, one speaks of `pure
Gaussian basis functions`:

.. math::
    P(r,\theta,\phi) = C_{\ell m}(r,\theta,\phi) \quad \text{or} \quad P(r,\theta,\phi) = S_{\ell m}(r,\theta,\phi)

where :math:`C_{\ell m}` and :math:`S_{\ell m}` are cosine and sine-like real regular
solid harmonics, defined as follows:

.. math::
    R_{\ell, 0}(r,\theta,\phi)  = C_{\ell 0}(r,\theta,\phi) & = R_\ell^0(r,\theta,\phi) \\
    R_{\ell, +m}(r,\theta,\phi) = C_{\ell m}(r,\theta,\phi) & = \frac{1}{\sqrt{2}}((-1)^m R_\ell^m(\theta,\phi) + R_\ell^{-m}(\theta,\phi)) \quad m = 1\ldots \ell \\
    R_{\ell, -m}(r,\theta,\phi) = S_{\ell m}(r,\theta,\phi) & = \frac{1}{i \sqrt{2}}((-1)^m R_\ell^m(\theta,\phi) - R_\ell^{-m}(\theta,\phi)) \quad m = 1\ldots \ell

where :math:`R_\ell^m` are the regular solid harmonics, which have in general
complex function values. The index :math:`\ell` is a zero or positive. Note that :math:`m` in
:math:`R_{\ell,+m}`, :math:`R_{\ell,-m}`, :math:`C_{\ell m}`, and :math:`S_{\ell
m}` is moved to the subscript to indicate that these are real functions. The
regular solid harmonics are derived from the standard spherical harmonics as
follows:

.. math::
    R_\ell^m(r, \theta, \varphi) & = \sqrt{\frac{4\pi}{2\ell+1}} r^l Y_\ell^m(\theta, \varphi) \\
        & = \sqrt{\frac{(\ell-m)!}{(\ell+m)!}} \, P_\ell^m(\cos{\theta})\, e^{i m \varphi}

(The Condon–Shortley phase is not used here.) Due to the presence of the factor
:math:`r^\ell`, the real regular solid harmonics are homogeneous polynomials,
i.e. linear combinations of the Cartesian polynomials defined in the previous
section.

Real regular solid harmonics are used because their normalization over the angular degrees of freedom, i.e.

.. math::
    \iint |R_{\ell m}(r, \theta, \varphi)|^2 \sin \theta d \theta d \varphi = \frac{4\pi r^{2l}}{2\ell+1},

is compatible with the Cartesian s- and p-type polynomials from the previous
section. The normalization constant of a pure Gaussian basis function is:

.. math:: N(\alpha_k, \ell) = \sqrt{\frac
        {(2\alpha_k/\pi)^{3/2} (4\alpha_k)^\ell}
        {(2\ell-1)!!}
        }

In practical applications, one always combines all the basis functions of a
given angular momentum. A basis specification typically only mentions the total
angular momentum, and it is assumed that all polynomials of that
order are included in the basis set. The number of basis functions, i.e. the
number of polynomials, for a given angular momentum, :math:`\ell`, is
:math:`2\ell+1`. For the implementation, one must fix a certain ordering of
these polynomials. The ordering in HORTON is based on the angular momentum
number, :math:`m`. When :math:`m>0`, the cosine-like functions is preceded by
the sine-like function.

The first five angular momenta and the corresponding polynomials are listed in
the table below.

====== ================ == ===========
Symbol Angular momentum #  Polynomials
====== ================ == ===========
S      0                1  :math:`C_{00}=1`
P      1                3  :math:`C_{10}=z`, :math:`C_{11}=x`, :math:`S_{11}=y`
D      2                5  :math:`C_{20}`, :math:`C_{21}`, :math:`S_{21}`, :math:`C_{22}`, :math:`S_{22}`
F      3                7  :math:`C_{30}`, :math:`C_{31}`, :math:`S_{31}`, :math:`C_{32}`, :math:`S_{32}`, :math:`C_{33}`, :math:`S_{33}`
G      4                9  :math:`C_{40}`, :math:`C_{41}`, :math:`S_{41}`, :math:`C_{42}`, :math:`S_{42}`, :math:`C_{43}`, :math:`S_{43}`, :math:`C_{44}`, :math:`S_{44}`
====== ================ == ===========


Transformation from Cartesian to pure basis functions
-----------------------------------------------------

Let us now derive convenient expressions for these real solid harmonics in terms
of Cartesian coordinates. The function :math:`P_\ell^m` is the
associated Legendre Polynomial. For positive :math:`m` we have:

.. math::
    P_\ell^m(x) & = (-1)^m (1-x^2)^{m/2} \frac{d^m}{dx^m} P_\ell(x) \\
    P_\ell^{-m}(x) & = (-1)^m \frac{(\ell-m)!}{(\ell+m)!} P_\ell^m

where :math:`P_\ell` is the ordinary Legendre polynomial of order :math:`\ell`.
Note that the factors :math:`(-1)^m` are canceled out in the definition of the
real solid harmonics. Substitution of these definitions leads to the following
form for the regular solid harmonics:

.. math::
    R_\ell^m(r, \theta, \varphi) = (-1)^{(m+|m|)/2}\sqrt{\frac{(\ell-|m|)!}{(\ell+|m|)!}} \, r^l sin^m \theta \frac{d^m P_\ell(\cos{\theta})}{d \cos \theta}\, e^{i m \varphi}

For :math:`m>0`, the real regular solid harmonics are first written as follows:

.. math::
    C_{\ell m}(r, \theta, \varphi) & = r^\ell \sqrt{\frac{2(\ell-m)!}{(\ell+m)!}} \, sin^m \theta \frac{d^m P_\ell(\cos{\theta})}{d \cos \theta}\, cos(m \varphi) \\
    S_{\ell m}(r, \theta, \varphi) & = r^\ell \sqrt{\frac{2(\ell-m)!}{(\ell+m)!}} \, sin^m \theta \frac{d^m P_\ell(\cos{\theta})}{d \cos \theta}\, sin(m \varphi)

It is conventional to factor out the :math:`z`-dependent part (which also has
some pure :math:`r`-dependence). Making use of :math:`z=r\cos\theta`, one gets:

.. math::
    \Pi_{\ell m}(z,r^2) & = r^{\ell-m} \frac{d^m P_\ell (\cos\theta)}{d \cos\theta} \\
             & = \sum_{k=0}^{\lfloor (\ell-m)/2 \rfloor} \gamma_{\ell k}^{(m)} r^{2k} z^{\ell-2k-m}

with

.. math::
    \gamma_{\ell k}^{(m)} = \frac{(-1)^k}{2^\ell} \binom{\ell}{k}\binom{2\ell-2k}{\ell}\frac{(\ell-2k)!}{(\ell-2k-m)!}

For the :math:`(x,y)`-dependence one has to define following polynomials for the
cosine and sine-like functions, respectively:

.. math::
    A_m(x,y) & = \mathrm{Re}[(x+iy)^m] \\
             & = r^m \sin^m \theta \cos(m \varphi) \\
             & = \frac{1}{2}\biggl( (r \sin \theta e^{i \varphi})^m + (r \sin \theta e^{-i \varphi})^m \biggr) \\
             & = \frac{1}{2}\biggl( (x + iy)^m + (x - iy)^m \biggr) \\
             & = \sum_p \binom{m}{p} x^p y^{m-p} \cos \bigl( (m-p) \pi/2 \bigl)

.. math::
    B_m(x,y) & = \mathrm{Im}[(x+iy)^m] \\
             & = r^m \sin^m \theta \sin(m \varphi) \\
             & = \frac{1}{2}\biggl( (r \sin \theta e^{i \varphi})^m - (r \sin \theta e^{-i \varphi})^m \biggr) \\
             & = \frac{1}{2}\biggl( (x + iy)^m - (x - iy)^m \biggr) \\
             & = \sum_p \binom{m}{p} x^p y^{m-p} \sin \bigl( (m-p) \pi/2 \bigl)

where we made use of :math:`i^k+(-i)^k = e^{k\pi/2} + e^{-k\pi/2} = \cos(k\pi/2)`
and :math:`i^k-(-i)^k = e^{k\pi/2} - e^{-k\pi/2} = \sin(k\pi/2)`. Putting it
all together, we have:

.. math::
    C_{\ell m}(x, y, z) & = \sqrt{\frac{2(\ell-m)!}{(\ell+m)!}} \, \Pi_{\ell m}(z,r^2) \, A_m(x,y) \\
    S_{\ell m}(x, y, z) & = \sqrt{\frac{2(\ell-m)!}{(\ell+m)!}} \, \Pi_{\ell m}(z,r^2) \, B_m(x,y)

Also for the case :math:`m=0`, one has a similar form:

.. math::
    C_{\ell 0}(x, y, z) & = \Pi_{\ell 0}(z,r^2) \\

These expressions allow one to write the real solid harmonics in terms of a
homogeneous polynomial of Cartesian coordinates. The following table is
generated by the script ``tools/harmonics.py``, which uses Sympy for the
symbolic manipulations:

.. math::
    C_{00}(x,y,z) & = 1 \\
    C_{10}(x,y,z) & = z \\
    C_{11}(x,y,z) & = x \\
    S_{11}(x,y,z) & = y \\
    C_{20}(x,y,z) & = - \frac{1}{2} r^{2} + \frac{3}{2} z^{2} \\
    C_{21}(x,y,z) & = \sqrt{3} x z \\
    S_{21}(x,y,z) & = \sqrt{3} y z \\
    C_{22}(x,y,z) & = \frac{1}{2} \sqrt{3} \left(x^{2} - y^{2}\right) \\
    S_{22}(x,y,z) & = \sqrt{3} x y \\
    C_{30}(x,y,z) & = - \frac{3}{2} r^{2} z + \frac{5}{2} z^{3} \\
    C_{31}(x,y,z) & = \frac{1}{6} \sqrt{6} x \left(- \frac{3}{2} r^{2} + \frac{15}{2} z^{2}\right) \\
    S_{31}(x,y,z) & = \frac{1}{6} \sqrt{6} y \left(- \frac{3}{2} r^{2} + \frac{15}{2} z^{2}\right) \\
    C_{32}(x,y,z) & = \frac{1}{2} \sqrt{15} z \left(x^{2} - y^{2}\right) \\
    S_{32}(x,y,z) & = \sqrt{15} x y z \\
    C_{33}(x,y,z) & = \frac{1}{4} \sqrt{10} \left(x^{3} - 3 x y^{2}\right) \\
    S_{33}(x,y,z) & = \frac{1}{4} \sqrt{10} \left(3 x^{2} y - y^{3}\right) \\
    C_{40}(x,y,z) & = \frac{3}{8} r^{4} - \frac{15}{4} r^{2} z^{2} + \frac{35}{8} z^{4} \\
    C_{41}(x,y,z) & = \frac{1}{10} \sqrt{10} x \left(- \frac{15}{2} r^{2} z + \frac{35}{2} z^{3}\right) \\
    S_{41}(x,y,z) & = \frac{1}{10} \sqrt{10} y \left(- \frac{15}{2} r^{2} z + \frac{35}{2} z^{3}\right) \\
    C_{42}(x,y,z) & = \frac{1}{30} \sqrt{5} \left(- \frac{15}{2} r^{2} + \frac{105}{2} z^{2}\right) \left(x^{2} - y^{2}\right) \\
    S_{42}(x,y,z) & = \frac{1}{15} \sqrt{5} x y \left(- \frac{15}{2} r^{2} + \frac{105}{2} z^{2}\right) \\
    C_{43}(x,y,z) & = \frac{1}{4} \sqrt{70} z \left(x^{3} - 3 x y^{2}\right) \\
    S_{43}(x,y,z) & = \frac{1}{4} \sqrt{70} z \left(3 x^{2} y - y^{3}\right) \\
    C_{44}(x,y,z) & = \frac{1}{8} \sqrt{35} \left(x^{4} - 6 x^{2} y^{2} + y^{4}\right) \\
    S_{44}(x,y,z) & = \frac{1}{8} \sqrt{35} \left(4 x^{3} y - 4 x y^{3}\right)


Note that these functions are not normalized yet.
The formatting of the list above is not great because of the limitations of
Sympy's latex printer.

The script ``tools/harmonics.py`` also generates the transformation matrices
from Cartesian to pure basis functions. These do take into account the
normalization.

.. math::
    \left(\begin{array}{c}
    X(C_{00})
    \end{array}\right)
    &=
    \left(\begin{array}{c}
    1 \\
    \end{array}\right)
    \left(\begin{array}{c}
    X(1)
    \end{array}\right)
    \\
    \left(\begin{array}{c}
    X(C_{10}) \\ X(C_{11}) \\ X(S_{11})
    \end{array}\right)
    &=
    \left(\begin{array}{ccc}
    0 & 0 & 1 \\
    1 & 0 & 0 \\
    0 & 1 & 0 \\
    \end{array}\right)
    \left(\begin{array}{c}
    X(x) \\ X(y) \\ X(z)
    \end{array}\right)
    \\
    \left(\begin{array}{c}
    X(C_{20}) \\ X(C_{21}) \\ X(S_{21}) \\ X(C_{22}) \\ X(S_{22})
    \end{array}\right)
    &=
    \left(\begin{array}{cccccc}
    - \frac{1}{2} & 0 & 0 & - \frac{1}{2} & 0 & 1 \\
    0 & 0 & 1 & 0 & 0 & 0 \\
    0 & 0 & 0 & 0 & 1 & 0 \\
    \frac{1}{2} \sqrt{3} & 0 & 0 & - \frac{1}{2} \sqrt{3} & 0 & 0 \\
    0 & 1 & 0 & 0 & 0 & 0 \\
    \end{array}\right)
    \left(\begin{array}{c}
    X(xx) \\ X(xy) \\ X(xz) \\ X(yy) \\ X(yz) \\ X(zz)
    \end{array}\right)
    \\
    \left(\begin{array}{c}
    X(C_{30}) \\ X(C_{31}) \\ X(S_{31}) \\ X(C_{32}) \\ X(S_{32}) \\ X(C_{33}) \\ X(S_{33})
    \end{array}\right)
    &=
    \left(\begin{array}{cccccccccc}
    0 & 0 & - \frac{3}{10} \sqrt{5} & 0 & 0 & 0 & 0 & - \frac{3}{10} \sqrt{5} & 0 & 1 \\
    - \frac{1}{4} \sqrt{6} & 0 & 0 & - \frac{1}{20} \sqrt{30} & 0 & \frac{1}{5} \sqrt{30} & 0 & 0 & 0 & 0 \\
    0 & - \frac{1}{20} \sqrt{30} & 0 & 0 & 0 & 0 & - \frac{1}{4} \sqrt{6} & 0 & \frac{1}{5} \sqrt{30} & 0 \\
    0 & 0 & \frac{1}{2} \sqrt{3} & 0 & 0 & 0 & 0 & - \frac{1}{2} \sqrt{3} & 0 & 0 \\
    0 & 0 & 0 & 0 & 1 & 0 & 0 & 0 & 0 & 0 \\
    \frac{1}{4} \sqrt{10} & 0 & 0 & - \frac{3}{4} \sqrt{2} & 0 & 0 & 0 & 0 & 0 & 0 \\
    0 & \frac{3}{4} \sqrt{2} & 0 & 0 & 0 & 0 & - \frac{1}{4} \sqrt{10} & 0 & 0 & 0 \\
    \end{array}\right)
    \left(\begin{array}{c}
    X(xxx) \\ X(xxy) \\ X(xxz) \\ X(xyy) \\ X(xyz) \\ X(xzz) \\ X(yyy) \\ X(yyz) \\ X(yzz) \\ X(zzz)
    \end{array}\right)
    \\
    \left(\begin{array}{c}
    X(C_{40}) \\ X(C_{41}) \\ X(S_{41}) \\ X(C_{42}) \\ X(S_{42}) \\ X(C_{43}) \\ X(S_{43}) \\ X(C_{44}) \\ X(S_{44})
    \end{array}\right)
    &=
    \left(\begin{array}{ccccccccccccccc}
    \frac{3}{8} & 0 & 0 & \frac{3}{140} \sqrt{105} & 0 & - \frac{3}{35} \sqrt{105} & 0 & 0 & 0 & 0 & \frac{3}{8} & 0 & - \frac{3}{35} \sqrt{105} & 0 & 1 \\
    0 & 0 & - \frac{3}{28} \sqrt{70} & 0 & 0 & 0 & 0 & - \frac{3}{28} \sqrt{14} & 0 & \frac{1}{7} \sqrt{70} & 0 & 0 & 0 & 0 & 0 \\
    0 & 0 & 0 & 0 & - \frac{3}{28} \sqrt{14} & 0 & 0 & 0 & 0 & 0 & 0 & - \frac{3}{28} \sqrt{70} & 0 & \frac{1}{7} \sqrt{70} & 0 \\
    - \frac{1}{4} \sqrt{5} & 0 & 0 & 0 & 0 & \frac{3}{14} \sqrt{21} & 0 & 0 & 0 & 0 & \frac{1}{4} \sqrt{5} & 0 & - \frac{3}{14} \sqrt{21} & 0 & 0 \\
    0 & - \frac{1}{14} \sqrt{35} & 0 & 0 & 0 & 0 & - \frac{1}{14} \sqrt{35} & 0 & \frac{3}{7} \sqrt{7} & 0 & 0 & 0 & 0 & 0 & 0 \\
    0 & 0 & \frac{1}{4} \sqrt{10} & 0 & 0 & 0 & 0 & - \frac{3}{4} \sqrt{2} & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
    0 & 0 & 0 & 0 & \frac{3}{4} \sqrt{2} & 0 & 0 & 0 & 0 & 0 & 0 & - \frac{1}{4} \sqrt{10} & 0 & 0 & 0 \\
    \frac{1}{8} \sqrt{35} & 0 & 0 & - \frac{3}{4} \sqrt{3} & 0 & 0 & 0 & 0 & 0 & 0 & \frac{1}{8} \sqrt{35} & 0 & 0 & 0 & 0 \\
    0 & \frac{1}{2} \sqrt{5} & 0 & 0 & 0 & 0 & - \frac{1}{2} \sqrt{5} & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
    \end{array}\right)
    \left(\begin{array}{c}
    X(xxxx) \\ X(xxxy) \\ X(xxxz) \\ X(xxyy) \\ X(xxyz) \\ X(xxzz) \\ X(xyyy) \\ X(xyyz) \\ X(xyzz) \\ X(xzzz) \\ X(yyyy) \\ X(yyyz) \\ X(yyzz) \\ X(yzzz) \\ X(zzzz)
    \end{array}\right)


These transformations are implemented in ``horton/cartpure.c`` with sparse
matrix products for angular momenta up to :math:`\ell=9`.


Recursion relations for real regular solid harmonics
----------------------------------------------------

Recurrence relations for :math:`\Pi_{\ell m}(z,r^2)` can be derived from the
recurrence relations for the associated Legendre polynomials:

    Initialization

    .. math::
        \Pi_{0 0}(z,r^2) & = 1

    For :math:`\ell \ge 1`

    .. math::
        \Pi_{\ell \ell}(z,r^2) & = (2\ell-1)\Pi_{\ell-1, \ell-1}(z,r^2) \\
        \Pi_{\ell, \ell-1}(z,r^2) & = z\Pi_{\ell\ell}(z,r^2)

    For :math:`\ell \ge 2` and :math:`0 \le m \le \ell-2`

    .. math::
        \Pi_{\ell m}(z,r^2) = z \frac{2\ell-1}{\ell-m} \Pi_{\ell-1, m}(z,r^2)
                            - r^2 \frac{\ell+m-1}{\ell-m} \Pi_{\ell-2, m}(z,r^2)


Recurrence relations for the functions :math:`A_m(x,y)` and :math:`B_m(x,y)` are
easily derived from scratch:

.. math::
    A_m(x,y) + i B_m(x,y) & = (x + iy)^m \\
                          & = (x + iy) (x + iy)^{m-1}\\
                          & = (x + iy) (A_{m-1}(x,y) + iB_{m-1}(x,y))

Hence, one gets:

    Initialization

        .. math::
            A_1(x,y) & = x \\
            B_1(x,y) & = y

    For :math:`m \ge 2`

        .. math::
            A_m(x,y) & = x A_{m-1}(x,y) - y B_{m-1}(x,y) \\
            B_m(x,y) & = y A_{m-1}(x,y) + x B_{m-1}(x,y)
