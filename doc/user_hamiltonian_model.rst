..
    : Horton is a development platform for electronic structure methods.
    : Copyright (C) 2011-2015 The Horton Development Team
    :
    : This file is part of Horton.
    :
    : Horton is free software; you can redistribute it and/or
    : modify it under the terms of the GNU General Public License
    : as published by the Free Software Foundation; either version 3
    : of the License, or (at your option) any later version.
    :
    : Horton is distributed in the hope that it will be useful,
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

Horton provides the following physical model Hamiltonians

1. The :ref:`1- dimensional Hubbard model Hamiltonian <hubbardham>` with open and periodic boundary conditions

Note that only ``DenseLinalgFactory`` is supported for model Hamiltonians.


.. _hubbardham:

The Hubbard model Hamiltonian
=============================

The Hubbard model Hamiltonian is the simplest model of interacting particles on a lattice and reads

.. math::
    :label: hubbard

    \hat{H}_{\rm Hub} = -t\sum_{j,\sigma} \left( a_{(j+1)\sigma}^{\dagger}a_{j\sigma}
    + a_{j\sigma}^{\dagger}a_{(j+1)\sigma} \right )
    +U\sum_j n_{j\uparrow} n_{j\downarrow},

where the first term is a one-electron term and accounts for nearest-neighbor hopping, while the second term is the repulsive on-site interaction. :math:`t` and :math:`U` are user specified parameters and :math:`\sigma` is the electron spin.


Preliminaries
-------------

In contrast to the molecular Hamiltonian, the Hubbard Hamiltonian does not require specifying a molecule and ``Gobasis`` instance. Instead, a ``LinalgFactory`` instance can be immediately created,

.. code-block:: python

   lf = DenseLinalgFactory(n)

where

   :lf: A linear algebra factory. Must be of type ``DenseLinalgFactory``. Note that ``CholeskyLinalgFactory`` is not supported
   :n: (int) number of lattice sites/orbitals


Defining the Hamiltonian and boundary conditions (open or periodic)
-------------------------------------------------------------------

Before the hopping and on-site repulsion terms can be calculated, an instance of the Hubbard Hamiltonian class ``Hubbard`` has to be created. To do so, run

.. code-block:: python

   modelham = Hubbard()

By default, the Hubbard Hamiltonian is constructed for periodic boundary conditions. Open boundary conditions can be enforced during initialization using

.. code-block:: python

   modelham = Hubbard(pbc=False)

with the only optional argument

    :pbc: (boolean) if ``True``, periodic boundary conditions are used (default ``True``)

The nearest-neighbor hopping term (:math:`t` in eq. :eq:`hubbard`) can be calculated as

.. code-block:: python

   hopping = modelham.compute_kinetic(lf, t)

where **lf** is again the linear algebra factory and

   :t: (float) the nearest-neighbor hopping parameter :math:`t`. Note that the hopping term is multiplied by a factor of -1 in the Hamiltonian (see eq. :eq:`hubbard`)

The strength of the on-site repulsion :math:`U` can be assigned as follows

.. code-block:: python

    onsite = modelham.compute_er(lf, U)

Similarly, **lf** is the linear algebra factory, and

   :U: (float) the on-site repulsion strength (:math:`U` parameter)

Finally, all terms of the 1-dimensinal Hubbard Hamiltonian can be combined together and passed to the effective Hamiltonian class ``REffHam``, which is used in the restricted Hartree-Fock model or DFT module,

.. code-block:: python

    terms = [
            RTwoIndexTerm(hopping, 'kin'),
            RDirectTerm(onsite, 'hartree'),
            RExchangeTerm(onsite, 'x_hf'),
            ]
    ham = REffHam(terms)

Note that the last step can be omitted for post-Hartree-Fock methods, like :ref:`AP1roG <introap1rog>`, :ref:`MP2 <mp2>`, :ref:`PTa <pta>`, and :ref:`PTb <ptb>`.


Filling the lattice
-------------------

The number of electrons/spin-1/2 particles on the lattice is assigned by the ``AufbauOccModel`` model. If the number of electrons/spin-1/2 particles is even and restricted orbitals are used, the occupation model can be set as follows

.. code-block:: python

   occ_model = AufbauOccModel(m)

where

   :m: (int) number of electron pairs/doubly occupied sites



The overlap matrix and initial guess orbitals
---------------------------------------------

To generate the initial orbitals for the `n` lattice sites, run

.. code-block:: python

    orb = lf.create_expansion(n)

where

   :n: (int) number of lattice sites/orbitals

The overlap matrix of the Hubbard model can be computed as follows

.. code-block:: python

    olp = modelham.compute_overlap(lf)

where **modelham** is an instance of ``Hubbard`` and **lf** is the linear algebra factory used.


Example Python script
=====================

Restricted Hartree-Fock calculations using the 1-dim Hubbard model Hamiltonian with PBC
---------------------------------------------------------------------------------------

This example shows a restricted Hartree-Fock calculation for the half-filled Hubbard model. Both the number of electron spins and sites is 6. The :math:`t` parameter is set to -1, while the :math:`U` parameter is equal to 2. Periodic boundary conditions are used.

.. literalinclude:: ../data/examples/hamiltonian/hubbard.py
    :caption: data/examples/hamiltonian/hubbard.py
    :lines: 2-
