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

.. _localization:

Localization of molecular orbitals
##################################

In general, the localization algorithm optimizes some localization function by
orthogonal transformation of the orbitals. Given orbitals
:math:`\vert i \rangle`, the localized orbital :math:`\vert \tilde{i} \rangle`
can be obtained by some transformation

.. math::

    \vert \tilde{i} \rangle = \sum_k \vert k \rangle \exp(-\mathbf{\kappa})_{ki},

where

.. math::

    \mathbf{\kappa} = \sum_{k > l} \kappa_{kl} (a^\dagger_k a_l - a^\dagger_l a_k)

and :math:`\kappa_{kl}` is determined by the optimization of the localization function.

Many of the localization schemes, and thus
the result of the localization, differ by the localization function. The localization
function somehow measures the localization of the orbitals. So far, HORTON only
supports the Pipek-Mezey localization. [pipek1989]_


.. _pipek-mezey:

Pipek-Mezey localization
========================

In the Pipek-Mezey scheme, the Pipek-Mezey localization function, :math:`D`, is maximized.

.. math::

    D = \sum_{i} \sum_{A \in \textrm{atoms}} (Q_{ii}^A)^2,

where :math:`\mathbf{Q}^A` is the atomic population matrix. The atomic population
matrix can be obtained from the overlap of the atomic basis, the expansion of the molecular
orbitals from the atomic basis, the occupation of each molecular orbital, and
some weighted projection of atomic basis function within each atom.

For example, if the Mulliken population analysis is used, the projectors are
obtained through :py:meth:`~horton.part.mulliken.get_mulliken_operators`.
Then the Pipek-Mezey localization function and the optimization are obtained through
:py:class:`~horton.localization.localization.PipekMezey`. Function call,
:py:meth:`~horton.localization.localization.Localization.__call__`, of this
instance results in localization.
Please see documentation in :py:mod:`horton.localization.localization` for more detail.


Example Python scripts
======================

Pipek-Mezey localization of restricted Hartree-Fock orbitals for the water molecule
-----------------------------------------------------------------------------------

This is a basic example on how to perform a Pipek-Mezey localization in HORTON. This script performs a Pipek-Mezey localization for the water molecule using the cc-pVDZ basis set and Mulliken projectors.


.. literalinclude :: ../data/examples/localization/water_pm.py
    :caption: data/examples/localization/water_pm.py
    :lines: 2-
