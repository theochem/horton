..
    : HORTON: Helpful Open-source Research TOol for N-fermion systems.
    : Copyright (C) 2011-2015 The HORTON Development Team
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

Orbital entanglement analysis
#############################

.. _orbital_entanglement:

Orbital entanglement and orbital correlation
============================================

In quantum chemistry, the interaction between orbitals is commonly used to
understand chemical processes. For example, molecular orbital diagrams, frontier
orbital theory, and ligand field theory all use orbitals to understand and justify
chemical phenomenon. The interaction between orbitals can be quantified through
ideas in entanglement.

Specifically, the interaction between one orbital and all the other orbitals can be
measured by the single-orbital entropy :math:`s(1)_i` (or one-orbital entropy).
It is calculated from the eigenvalues :math:`\omega_{\alpha,i}` of the
one-orbital reduced density matrix (RDM)

.. math::

    s(1)_i = -\sum_{\alpha=1}^4 \omega_{\alpha,i}\ln \omega_{\alpha,i}.

The one-orbital RDM is obtained from the N-particle RDM by tracing out all
other orbital degrees of freedom except those of orbital *i*. This leads to an
RDM whose dimension is equal to the dimension of the one-orbital Fock space. For
spatial orbitals, the one-orbital Fock space has a dimension of 4 (one for unoccupied,
one for doubly occupied and two for singly occupied).

Similarly, the two-orbital entropy :math:`s(2)_{i,j}` quantifies the interaction
of an orbital pair :math:`i,j` and all other orbitals. It is calculated from the
eigenvalues of the two-orbital RDM :math:`\omega_{\alpha, i, j}` (with 16
possible states for spatial orbitals)

.. math::

    s(2)_{i,j} =-\sum_{\alpha=1}^{16} \omega_{\alpha, i, j} \ln \omega_{\alpha, i, j}.

The total amount of correlation between any pair of orbitals :math:`(i,j)` can
be evaluated from the orbital-pair mutual information. The orbital-pair mutual
information is calculated using the single- and two-orbital entropy and thus
requires the one- and two-orbital RDMs,

.. math::

    I_{i|j} = \frac{1}{2} \big(s(2)_{i,j} - s(1)_{i} - s(1)_{j} \big) \big(1 - \delta_{ij}\big),

where :math:`\delta_{ij}` is the Kronecker delta.

Note that a correlated wavefunction is required to have non-zero orbital entanglement
and correlation. In the case of an uncorrelated wavefunction (for instance, a
single Slater determinant) the (orbital) entanglement entropy is zero.

For more information on orbital entanglement and correlation, see ref. [boguslawski2015a]_.


Supported features
==================

Unless mentioned otherwise, the orbital entanglement module only supports restricted
orbitals, ``DenseLinalgFactory`` and ``CholeskyLinalgFactory``. The current
version of HORTON offers

1. :ref:`Calculation of single-orbital entropy and orbital-pair mutual
information for a seniority-zero wavefunction <orbital_entanglementseniorityzero>`

Supported wavefunction models are

1. :ref:`AP1roG <user_ap1rog>` (seniority-zero wavefunction)


.. _orbital_entanglementseniorityzero:

Orbital entanglement and correlation for a seniority zero wavefunction
======================================================================

If you use this module, please cite [boguslawski2015a]_

How to set-up a calculation
---------------------------

To evaluate the single-orbital entropy and orbital-pair mutual information for a
given wavefunction model, the corresponding one and two-particle RDMs need to be calculated
first. Then one- and two-particle RDMs are used to evaluate the one- and two-orbital RDM's
from which the eigenvalues are needed to calculate the single-orbital and
two-orbital entropy.

Given the one- and two-particle RDMs, an instance of the
:py:class:`~horton.orbital_entanglement.orbital_entanglement.orbital_entanglement`
can be created. A function call of this instance,
:py:meth:`~horton.orbital_entanglement.orbital_entanglement.orbital_entanglement.__call__`,
calculates the entanglement and correlation measures.


How to generate correlation diagrams
------------------------------------

To generate the single-orbital entropy diagram and mutual information plot, run
in the terminal:

.. code-block:: bash

    horton-entanglement.py threshold init_index final_index

where **threshold** determines the lower cutoff value of the mutual information
and must be given in orders of magnitude (1, 0.1, 0.01, 0.001, etc.). Orbital
correlations that are smaller than **cutoff** will not be displayed in the
mutual information diagram. The arguments **init_index** and **final_index** are
optional arguments. If provided, the single-orbital entropy will be plotted for
orbital indices in the interval [init_index, final_index].


Example Python scripts
======================

Orbital entanglement analysis of an AP1roG wavefunction
-------------------------------------------------------

This is a basic example on how to perform an orbital entanglement analysis in
HORTON. This script performs an orbital-optimized AP1roG calculation, followed
by an orbital entanglement analysis of the AP1roG wavefunction for the water
molecule using the cc-pVDZ basis set.

.. literalinclude :: ../data/examples/orbital_entanglement/water.py
    :caption: data/examples/orbital_entanglement/water.py
    :lines: 2-
