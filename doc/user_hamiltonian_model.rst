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

.. _modphysham:

Defining a model Hamiltonian
############################

Supported features
==================

HORTON implements the :ref:`1- dimensional Hubbard model Hamiltonian
<hubbardham>` with open and periodic boundary conditions. Note that only
``DenseLinalgFactory`` is supported for model Hamiltonians.

.. _hubbardham:

The Hubbard model Hamiltonian
=============================

The Hubbard model Hamiltonian is the simplest model of interacting particles on a lattice and reads

.. math::
    :label: hubbard

    \hat{H}_{\rm Hub} = -t\sum_{j,\sigma} \left( a_{(j+1)\sigma}^{\dagger}a_{j\sigma}
    + a_{j\sigma}^{\dagger}a_{(j+1)\sigma} \right )
    +U\sum_j n_{j\uparrow} n_{j\downarrow},

where the first term is a one-electron term and accounts for the nearest-neighbor
hopping, while the second term is the repulsive on-site interaction. The :math:`t`
and :math:`U` are user specified parameters and :math:`\sigma` is the electron spin.


Preliminaries
-------------

In contrast to the molecular Hamiltonian, the Hubbard Hamiltonian does not require
a molecule and :py:class:`~horton.gbasis.cext.GOBasis` instance. Instead, a
:py:class:`~horton.matrix.dense.DenseLinalgFactory` instance can
be immediately created,

.. code-block:: python

   lf = DenseLinalgFactory(n)

Note that ``CholeskyLinalgFactory`` is not supported.


Defining the Hamiltonian and boundary conditions (open or periodic)
-------------------------------------------------------------------

Before the hopping and on-site repulsion terms can be calculated, you need to create
an instance of :py:class:`~horton.modelhamiltonians.physmodham.Hubbard`:

.. code-block:: python

   modelham = Hubbard()

By default, the Hubbard Hamiltonian is constructed with periodic boundary conditions.
Open boundary conditions can be enforced during the initialization using

.. code-block:: python

   modelham = Hubbard(pbc=False)

The nearest-neighbor hopping term (:math:`t` in eq. :eq:`hubbard`) is calculated
using the method :py:meth:`~horton.modelhamiltonians.physmodham.Hubbard.compute_kinetic`:

.. code-block:: python

   hopping = modelham.compute_kinetic(lf, t)

The the on-site repulsion :math:`U` is calculated using the method
:py:meth:`~horton.modelhamiltonians.physmodham.Hubbard.compute_er`:

.. code-block:: python

    onsite = modelham.compute_er(lf, U)

Finally, all terms of the 1-dimensional Hubbard Hamiltonian are combined together
and passed to the effective Hamiltonian class, :py:class:`~horton.meanfield.hamiltonian.REffHam`,
which can then be passed to the restricted Hartree-Fock or DFT modules,

.. code-block:: python

    terms = [
            RTwoIndexTerm(hopping, 'kin'),
            RDirectTerm(onsite, 'hartree'),
            RExchangeTerm(onsite, 'x_hf'),
            ]
    ham = REffHam(terms)


Filling the lattice
-------------------

The number of electrons/spin-half particles on the lattice is assigned using the
Aufbau occupation number model. An instance of the class
:py:class:`~horton.meanfield.occ.AufbauOccModel` will be used later for setting
the occupations.

.. code-block:: python

   occ_model = AufbauOccModel(m)

Note that the number of electrons/spin-half particles must be even and orbitals
must be restricted.


Generate initial guess orbitals and overlap matrix
--------------------------------------------------

The instance of class :py:class:`~horton.matrix.dense.DenseExpansion` by using
the method :py:meth:`~horton.matrix.dense.DenseLinalgFactory.create_expansion`
is used to initiate the orbitals:

.. code-block:: python

    orb = lf.create_expansion(n)

An initial guess can be obtained with
:py:func:`~horton.meanfield.guess.guess_core_hamiltonian`:

.. code-block:: python

    guess_core_hamiltonian(olp, hopping, orb)

The overlap matrix of the Hubbard model can be computed using the method
:py:meth:`~horton.modelhamiltonians.physmodham.Hubbard.compute_overlap`:

.. code-block:: python

    olp = modelham.compute_overlap(lf)


Example Python script
=====================

Restricted Hartree-Fock calculations using the 1-dim Hubbard model Hamiltonian with PBC
---------------------------------------------------------------------------------------

This example shows a restricted Hartree-Fock calculation for the half-filled
Hubbard model. Both the number of electron spins and sites is 6. The :math:`t`
parameter is set to -1, while the :math:`U` parameter is equal to 2. Periodic
boundary conditions are used.

.. literalinclude:: ../data/examples/hamiltonian/hubbard.py
    :caption: data/examples/hamiltonian/hubbard.py
    :lines: 3-38
