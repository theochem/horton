
Gaussian basis functions
########################


Introduction
============


Horton supports contracted Gaussian basis functions, which have in general the
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
that are needed for the implementation in Horton.


Cartesian
=========


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
these polynomials. In Horton, the ordering is simply alphabetical.

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


Pure or Harmonic
================

When the polynomial is a real regular solid harmonic, one speaks of `pure
Gaussian basis functions`:

.. math::
    P(r,\theta,\phi) = C_\ell^m(r,\theta,\phi) \quad \text{or} \quad P(r,\theta,\phi) = S_\ell^m(r,\theta,\phi)

where :math:`C_\ell^m` and :math:`S_\ell^m` are cosine and sine-like real regular
solid harmonics, defined as follows:

.. math::
    C_\ell^0(r,\theta,\phi) & = r^\ell Y_\ell^0(r,\theta,\phi) \\
    C_\ell^m(r,\theta,\phi) & = r^\ell \frac{1}{\sqrt{2}}((-1)^m Y_\ell^m(\theta,\phi) + Y_\ell^{-m}(\theta,\phi)) \quad m = 1\ldots \ell \\
    S_\ell^m(r,\theta,\phi) & = r^\ell \frac{1}{i \sqrt{2}}((-1)^m Y_\ell^m(\theta,\phi) - Y_\ell^{-m}(\theta,\phi)) \quad m = 1\ldots \ell

The index :math:`\ell` is a zero or positive. Due to the presence of the factor
:math:`r^\ell`, these are homogeneous polynomials, i.e. linear combinations of the
Cartesian polynomials listed above. The following definition for the spherical
harmonics is used:

.. math::
    Y_\ell^m(\theta, \varphi) = \sqrt{\frac{(\ell-m)!}{(\ell+m)!}} \, P_\ell^m(\cos{\theta})\, e^{i m \varphi}

which has a norm

.. math::
    \iint |Y_\ell^m(\theta, \varphi)|^2 \sin \theta d \theta d \varphi = \frac{4\pi}{2\ell+1}

for the sake of compatibility with the unnormalized Cartesian s- and p-type
functions above. The normalization constant of a pure Gaussian basis function
is:

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
these polynomials. The ordering in Horton is based on the angular momentum
number, :math:`m`. When :math:`m>0`, the cosine-like functions is preceded by
the sine-like function.

The first five angular momenta and the corresponding polynomials are listed in
the table below.

====== ================ == ===========
Symbol Angular momentum #  Polynomials
====== ================ == ===========
S      0                1  :math:`C_0^0=1`
P      1                3  :math:`C_1^0=z`, :math:`C_1^1=x`, :math:`S_1^1=y`
D      2                8  :math:`C_2^0`, :math:`C_2^1`, :math:`S_2^1`, :math:`C_2^2`, :math:`S_2^2`
F      3                7  :math:`C_3^0`, :math:`C_3^1`, :math:`S_3^1`, :math:`C_3^2`, :math:`S_3^2`, :math:`C_3^3`, :math:`S_3^3`
G      4                9  :math:`C_4^0`, :math:`C_4^1`, :math:`S_4^1`, :math:`C_4^2`, :math:`S_4^2`, :math:`C_4^3`, :math:`S_4^3`, :math:`C_4^4`, :math:`S_4^4`
====== ================ == ===========


Transformation from Cartesian to pure
=====================================

Let us now derive convenient expressions for these real solid harmonics in terms
of Cartesian coordinates. The function :math:`P_\ell^m` is the
associated Legendre Polynomial. For positive :math:`m` we have:

.. math::
    P_\ell^m(x) & = (-1)^m (1-x^2)^{m/2} \frac{d^m}{dx^m} P_\ell(x) \\
    P_\ell^{-m}(x) & = (-1)^m \frac{(\ell-m)!}{(\ell+m)!} P_\ell^m

where :math:`P_\ell` is the ordinary Legendre polynomial of order :math:`\ell`
Note that the factors :math:`(-1)^m` are compensated in the definition of the
real solid harmonics. Substitution of these definition leads to the following
form for the spherical harmonics:

.. math::
    Y_\ell^m(\theta, \varphi) = (-1)^{(m+|m|)/2}\sqrt{\frac{(\ell-|m|)!}{(\ell+|m|)!}} \, sin^m \theta \frac{d^m P_\ell(\cos{\theta})}{d \cos \theta}\, e^{i m \varphi}

For :math:`m=0`, it is trivial to write the Cartesian form of the real solid
harmonic with the substitution :math:`z=r\cos\theta`:

.. math::
    C_\ell^0(x, y, z) = P_\ell(z)

For :math:`m>0`, the real spherical harmonics are first written as follows:

.. math::
    C_\ell^m(\theta, \varphi) & = r^\ell \sqrt{\frac{2(\ell-m)!}{(\ell+m)!}} \, sin^m \theta \frac{d^m P_\ell(\cos{\theta})}{d \cos \theta}\, cos(m \varphi) \\
    S_\ell^m(\theta, \varphi) & = r^\ell \sqrt{\frac{2(\ell-m)!}{(\ell+m)!}} \, sin^m \theta \frac{d^m P_\ell(\cos{\theta})}{d \cos \theta}\, sin(m \varphi)

It is conventional to factor out the :math:`z`-dependent part (which also has
some pure :math:`r`-dependence). Making use of :math:`z=r\cos\theta`, one gets:

.. math::
    \Pi_\ell^m(z,r^2) & = r^{\ell-m} \frac{d^m P_\ell (\cos\theta)}{d \cos\theta} \\
             & = \sum_{k=0}^{\lfloor (\ell-m)/2 \rfloor} \gamma_{\ell k}^{(m)} r^{2k} z^{\ell-2k-m}

with

.. math::
    \gamma_{\ell k}^{(m)} = \frac{(-1)^k}{2^\ell} \binom{\ell}{k}\binom{2\ell-2k}{\ell}\frac{(\ell-2k)!}{(\ell-2k-m)!}

For the :math:`(x,y)`-dependence one has to following polynomials for the
cosine and sine-like functions, respectively:

.. math::
    A_m(x,y) & = r^m \sin^m \theta \cos(m \varphi) \\
             & = \frac{1}{2}\biggl( (r \sin \theta e^{i \varphi})^m + (r \sin \theta e^{-i \varphi})^m \biggr) \\
             & = \frac{1}{2}\biggl( (x + iy)^m + (x - iy)^m \biggr) \\
             & = \sum_p \binom{m}{p} x^p y^m-p \cos \bigl( (m-p) \pi/2 \bigl)

.. math::
    B_m(x,y) & = r^m \sin^m \theta \sin(m \varphi) \\
             & = \frac{1}{2}\biggl( (r \sin \theta e^{i \varphi})^m - (r \sin \theta e^{-i \varphi})^m \biggr) \\
             & = \frac{1}{2}\biggl( (x + iy)^m - (x - iy)^m \biggr) \\
             & = \sum_p \binom{m}{p} x^p y^m-p \sin \bigl( (m-p) \pi/2 \bigl)

where we made use of :math:`i^k+(-i)^k = e^{k\pi/2} + e^{-k\pi/2} = \cos(k\pi/2)`
and :math:`i^k-(-i)^k = e^{k\pi/2} - e^{-k\pi/2} = \sin(k\pi/2)`. Putting it
all together, we have:

.. math::
    C_\ell^m(x, y, z) & = \sqrt{\frac{2(\ell-m)!}{(\ell+m)!}} \, \Pi_\ell^m(z,r^2) \, A_m(x,y) \\
    S_\ell^m(x, y, z) & = \sqrt{\frac{2(\ell-m)!}{(\ell+m)!}} \, \Pi_\ell^m(z,r^2) \, B_m(x,y)

Also for the case :math:`m=0`, one has a similar form:

.. math::
    C_\ell^0(x, y, z) & = \Pi_\ell^0(z,r^2) \\

These expressions allow one to write the real solid harmonics in terms of a
homogeneous polynomial of Cartesian coordinates. The following table is
generated by the script ``tools/harmonics.py``, which uses Sympy for the
symbolic manipulations:

.. math::
    C_0^0(x,y,z) & = 1 \\
    C_1^0(x,y,z) & = z \\
    C_1^1(x,y,z) & = x \\
    S_1^1(x,y,z) & = y \\
    C_2^0(x,y,z) & = - \frac{1}{2} r^{2} + \frac{3}{2} z^{2} \\
    C_2^1(x,y,z) & = \sqrt{3} x z \\
    S_2^1(x,y,z) & = \sqrt{3} y z \\
    C_2^2(x,y,z) & = \frac{1}{2} \sqrt{3} \left(x^{2} - y^{2}\right) \\
    S_2^2(x,y,z) & = \sqrt{3} x y \\
    C_3^0(x,y,z) & = - \frac{3}{2} r^{2} z + \frac{5}{2} z^{3} \\
    C_3^1(x,y,z) & = \frac{1}{6} \sqrt{6} x \left(- \frac{3}{2} r^{2} + \frac{15}{2} z^{2}\right) \\
    S_3^1(x,y,z) & = \frac{1}{6} \sqrt{6} y \left(- \frac{3}{2} r^{2} + \frac{15}{2} z^{2}\right) \\
    C_3^2(x,y,z) & = \frac{1}{2} \sqrt{15} z \left(x^{2} - y^{2}\right) \\
    S_3^2(x,y,z) & = \sqrt{15} x y z \\
    C_3^3(x,y,z) & = \frac{1}{4} \sqrt{10} \left(x^{3} - 3 x y^{2}\right) \\
    S_3^3(x,y,z) & = \frac{1}{4} \sqrt{10} \left(3 x^{2} y - y^{3}\right) \\
    C_4^0(x,y,z) & = \frac{3}{8} r^{4} - \frac{15}{4} r^{2} z^{2} + \frac{35}{8} z^{4} \\
    C_4^1(x,y,z) & = \frac{1}{10} \sqrt{10} x \left(- \frac{15}{2} r^{2} z + \frac{35}{2} z^{3}\right) \\
    S_4^1(x,y,z) & = \frac{1}{10} \sqrt{10} y \left(- \frac{15}{2} r^{2} z + \frac{35}{2} z^{3}\right) \\
    C_4^2(x,y,z) & = \frac{1}{30} \sqrt{5} \left(- \frac{15}{2} r^{2} + \frac{105}{2} z^{2}\right) \left(x^{2} - y^{2}\right) \\
    S_4^2(x,y,z) & = \frac{1}{15} \sqrt{5} x y \left(- \frac{15}{2} r^{2} + \frac{105}{2} z^{2}\right) \\
    C_4^3(x,y,z) & = \frac{1}{4} \sqrt{70} z \left(x^{3} - 3 x y^{2}\right) \\
    S_4^3(x,y,z) & = \frac{1}{4} \sqrt{70} z \left(3 x^{2} y - y^{3}\right) \\
    C_4^4(x,y,z) & = \frac{1}{8} \sqrt{35} \left(x^{4} - 6 x^{2} y^{2} + y^{4}\right) \\
    S_4^4(x,y,z) & = \frac{1}{8} \sqrt{35} \left(4 x^{3} y - 4 x y^{3}\right)


Note that these functions are not normalized yet.
The formatting of the list above is not great because of the limitations of
Sympy's latex printer.

The script ``tools/harmonics.py`` also generates the transformation matrices
from Cartesian to pure basis functions. These do take into account the
normalization.

.. math::
    \left(\begin{array}{c}
    X(C_0^0)
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
    X(C_1^0) \\ X(C_1^1) \\ X(S_1^1)
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
    X(C_2^0) \\ X(C_2^1) \\ X(S_2^1) \\ X(C_2^2) \\ X(S_2^2)
    \end{array}\right)
    &=
    \left(\begin{array}{cccccc}
    -1 & 0 & 0 & - \frac{1}{2} & 0 & \frac{3}{2} \\
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
    X(C_3^0) \\ X(C_3^1) \\ X(S_3^1) \\ X(C_3^2) \\ X(S_3^2) \\ X(C_3^3) \\ X(S_3^3)
    \end{array}\right)
    &=
    \left(\begin{array}{cccccccccc}
    0 & 0 & - \frac{3}{5} \sqrt{5} & 0 & 0 & 0 & 0 & - \frac{3}{10} \sqrt{5} & 0 & \frac{5}{2} \\
    - \frac{1}{2} \sqrt{6} & 0 & 0 & - \frac{1}{20} \sqrt{30} & 0 & \frac{1}{4} \sqrt{30} & 0 & 0 & 0 & 0 \\
    0 & - \frac{1}{10} \sqrt{30} & 0 & 0 & 0 & 0 & - \frac{1}{4} \sqrt{6} & 0 & \frac{1}{4} \sqrt{30} & 0 \\
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
    X(C_4^0) \\ X(C_4^1) \\ X(S_4^1) \\ X(C_4^2) \\ X(S_4^2) \\ X(C_4^3) \\ X(S_4^3) \\ X(C_4^4) \\ X(S_4^4)
    \end{array}\right)
    &=
    \left(\begin{array}{ccccccccccccccc}
    \frac{3}{2} & 0 & 0 & \frac{3}{70} \sqrt{105} & 0 & - \frac{3}{14} \sqrt{105} & 0 & 0 & 0 & 0 & \frac{3}{8} & 0 & - \frac{3}{28} \sqrt{105} & 0 & \frac{35}{8} \\
    0 & 0 & - \frac{3}{14} \sqrt{70} & 0 & 0 & 0 & 0 & - \frac{3}{28} \sqrt{14} & 0 & \frac{1}{4} \sqrt{70} & 0 & 0 & 0 & 0 & 0 \\
    0 & 0 & 0 & 0 & - \frac{3}{14} \sqrt{14} & 0 & 0 & 0 & 0 & 0 & 0 & - \frac{3}{28} \sqrt{70} & 0 & \frac{1}{4} \sqrt{70} & 0 \\
    - \frac{1}{2} \sqrt{5} & 0 & 0 & \frac{1}{28} \sqrt{21} & 0 & \frac{1}{4} \sqrt{21} & 0 & 0 & 0 & 0 & \frac{1}{4} \sqrt{5} & 0 & - \frac{1}{4} \sqrt{21} & 0 & 0 \\
    0 & - \frac{1}{7} \sqrt{35} & 0 & 0 & 0 & 0 & - \frac{1}{14} \sqrt{35} & 0 & \frac{1}{2} \sqrt{7} & 0 & 0 & 0 & 0 & 0 & 0 \\
    0 & 0 & \frac{1}{4} \sqrt{10} & 0 & 0 & 0 & 0 & - \frac{3}{4} \sqrt{2} & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
    0 & 0 & 0 & 0 & \frac{3}{4} \sqrt{2} & 0 & 0 & 0 & 0 & 0 & 0 & - \frac{1}{4} \sqrt{10} & 0 & 0 & 0 \\
    \frac{1}{8} \sqrt{35} & 0 & 0 & - \frac{3}{4} \sqrt{3} & 0 & 0 & 0 & 0 & 0 & 0 & \frac{1}{8} \sqrt{35} & 0 & 0 & 0 & 0 \\
    0 & \frac{1}{2} \sqrt{5} & 0 & 0 & 0 & 0 & - \frac{1}{2} \sqrt{5} & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
    \end{array}\right)
    \left(\begin{array}{c}
    X(xxxx) \\ X(xxxy) \\ X(xxxz) \\ X(xxyy) \\ X(xxyz) \\ X(xxzz) \\ X(xyyy) \\ X(xyyz) \\ X(xyzz) \\ X(xzzz) \\ X(yyyy) \\ X(yyyz) \\ X(yyzz) \\ X(yzzz) \\ X(zzzz)
    \end{array}\right)
    \\
