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

The perturbation theory module
##############################

Supported features
==================

The perturbation theory module supports spin-restricted orbitals and the ``DenseLinalgFactory``. Horton offers the following flavors of perturbation theory:

1. Moller-Plesset Perturbation theory of second order with a restricted, closed-shell Hartree-Fock reference function (see :ref:`mp2`)

2. PTa with an AP1roG reference function (see :ref:`pta`)

3. PTb with an AP1roG reference function (see :ref:`ptb`)


.. _mp2:

Moller-Plesset perturbation theory of second order
==================================================

.. _getstartedmp2:

Getting started
---------------

The RMP2 module requires a restricted Hartree-Fock reference wavefunction (see :ref:`user_hf_dft`), its energy, and the one- and two-electron integrals as input arguments.

The Hartree-Fock energy can be calculated after SCF convergence as follows

.. code-block:: python

    ehf = ham.compute_energy()

where ``ham`` is an instance of the restricted effective Hamiltonian class ``REffHam`` (see FIXME), which contains all one- and two-electron terms. Note that all one-electron terms have to be combined into one single operator term. This can be done in the following way if the Hamiltonian contains the kinetic energy of the electrons and the electron-nuclear attraction,

.. code-block:: python

    ###############################################################################
    ## Calculate kinetic energy (kin) and nuclear attraction (na) term ############
    ###############################################################################
    kin = obasis.compute_kinetic(lf)
    na = obasis.compute_nuclear_attraction(mol.coordinates, mol.pseudo_numbers, lf)
    ###############################################################################
    ## Combine one-electron integrals to single Hamiltonian #######################
    ###############################################################################
    one = kin.copy()
    one.iadd(na)


How to set-up an MP2 calculation
--------------------------------

First, create an instance of the ``RMP2`` class

.. code-block:: python

    mp2 = RMP2(lf, occ_model)

with arguments

    :lf: A ``LinalgFactory`` instance (see FIXME). Only ``DenseLinalgFactory`` is supported
    :occ_model: (``AufbauOccModel`` instance) an Aufbau occupation model

A function call initiates an MP2 calculation,

.. code-block:: python

    emp2, tmp2 = mp2(one, two, orb, **{'eref': ehf, 'indextrans': 'tensordot'})

with arguments

    :one: (``TwoIndex`` instance) the one-electron integrals
    :two: (``FourIndex`` instance) the two-electron integrals
    :orb: (``Expansion`` instance) the AO/MO coefficient matrix

and keyword arguments

    :eref: (float) the Hartree-Fock reference energy (default ``float('nan')``) (see :ref:`getstartedmp2`)
    :indextrans: (str, optional) the 4-index transformation. Choice between ``tensordot`` (default) and ``einsum``. ``tensordot`` is faster than ``einsum``, requires, however, more memory. If ``DenseLinalgFactory`` is used, the memory requirement scales as :math:`2N^4` for ``einsum`` and :math:`3N^4` for ``tensordot``, respectively. Due to the storage of the two-electron integrals, the total amount of memory increases to :math:`3N^4` for ``einsum`` and :math:`4N^4` for ``tensordot``, respectively.

The function call gives 2 return values:

    :emp2: (list of float) the MP2 energy correction (first element)
    :tmp2: (list of ``FourIndex`` instances) the MP2 amplitudes (first element). The double excitation amplitudes :math:`t_{ij}^{ab}` are stored as ``t[i,a,j,b]``


.. _pta:

The PTa module
==============

.. _getstartedpta:

Getting started
---------------

PTa [Limacher2014]_ adds dynamic electron correlation effects on top of an AP1roG wavefunction (see :ref:`introap1rog`). You have to optimize an AP1roG wavefunction (see :ref:`ooap1rog` or :ref:`ap1rog`), before the PTa energy correction can be determined.


How to set-up a calculation
---------------------------

First, create an instance of the ``PTa`` class

.. code-block:: python

    pta = PTa(lf, occ_model)

with arguments

    :lf: A ``LinalgFactory`` instance (see FIXME). Only ``DenseLinalgFactory`` is supported
    :occ_model: (``AufbauOccModel`` instance) an Aufbau occupation model

A function call initiates an PTa calculation,

.. code-block:: python

    epta, tpta = pta(one, two, orb, c, **{'eref': energy, 'ecore': ecore, 'indextrans': 'tensordot'})

with arguments

    :one: (``TwoIndex`` instance) the one-electron integrals (the same integrals as used in the AP1roG module :ref:`introap1rog`)
    :two: (``FourIndex`` instance) the two-electron integrals (the same integrals as used in the AP1roG module :ref:`introap1rog`)
    :orb: (``Expansion`` instance) the optimized AP1roG MO coefficient matrix
    :c: (``TwoIndex`` instance) the geminal coefficients (see :ref:`ooap1rog`)

and keyword arguments

    :eref: (float) the AP1roG reference energy (default ``float('nan')``) (see :ref:`ooap1rog` how to get the AP1roG reference energy)
    :ecore: (float) the core energy (default ``float('nan')``). Usually, the nuclear repulsion term
    :indextrans: (str, optional) the 4-index transformation. Choice between ``tensordot`` (default) and ``einsum``. ``tensordot`` is faster than ``einsum``, requires, however, more memory. If ``DenseLinalgFactory`` is used, the memory requirement scales as :math:`2N^4` for ``einsum`` and :math:`3N^4` for ``tensordot``, respectively. Due to the storage of the two-electron integrals, the total amount of memory increases to :math:`3N^4` for ``einsum`` and :math:`4N^4` for ``tensordot``, respectively.

The function call gives 2 return values:

    :epta: (list of float) the PTa energy corrections. Contains the total PTa energy correction (first element), its seniority-zero contribution (second element), its seniority-two contribution (third element), and its seniority-four contribution (fourth element)
    :tpta: (list of ``FourIndex`` instances) the PTa amplitudes (first element). The double excitation amplitudes :math:`t_{ij}^{ab}` are stored as ``t[i,a,j,b]``


.. _ptb:

The PTb module
==============

.. _getstartedptb:

Getting started
---------------

PTb [Limacher2014]_ represents a different flavor to add dynamic electron correlation effects on top of an AP1roG wavefunction (see :ref:`introap1rog`). Similarly to PTa (:ref:`pta`), you have to optimize an AP1roG wavefunction (see :ref:`ooap1rog` or :ref:`ap1rog`), before the PTb energy correction can be determined.

How to set-up a calculation
---------------------------

First, create an instance of the ``PTb`` class

.. code-block:: python

    ptb = PTb(lf, occ_model)

with arguments

    :lf: A ``LinalgFactory`` instance (see FIXME). Only ``DenseLinalgFactory`` is supported
    :occ_model: (``AufbauOccModel`` instance) an Aufbau occupation model

A function call initiates an PTb calculation,

.. code-block:: python

    eptb, tptb = ptb(one, two, orb, c, **{'eref': energy, 'ecore': ecore})

with arguments

    :one: (``TwoIndex`` instance) the one-electron integrals (the same integrals as used in the AP1roG module :ref:`introap1rog`)
    :two: (``FourIndex`` instance) the two-electron integrals (the same integrals as used in the AP1roG module :ref:`introap1rog`)
    :orb: (``Expansion`` instance) the optimized AP1roG MO coefficient matrix
    :c: (``TwoIndex`` instance) the geminal coefficients (see :ref:`ooap1rog`)

Note that optional keyword arguments have been omitted. All keyword arguments are summarized in :ref:`ptbkeywords`.

The function call gives 2 return values:

    :eptb: (list of float) the PTb energy corrections. Contains the total PTb energy correction (first element), its seniority-zero contribution (second element), its seniority-two contribution (third element), and its seniority-four contribution (fourth element)
    :tptb: (list of ``FourIndex`` instances) the PTb amplitudes (first element). The double excitation amplitudes :math:`t_{ij}^{ab}` are stored as ``t[i,a,j,b]``

.. _ptbkeywords:


Summary of keyword arguments
----------------------------

    :indextrans: (str) 4-index Transformation. Choice between ``tensordot`` (default) and ``einsum``. ``tensordot`` is faster than ``einsum``, requires, however, more memory. If ``DenseLinalgFactory`` is used, the memory requirement scales as :math:`2N^4` for ``einsum`` and :math:`3N^4` for ``tensordot``, respectively. Due to the storage of the two-electron integrals, the total amount of memory increases to :math:`3N^4` for ``einsum`` and :math:`4N^4` for ``tensordot``, respectively.

    :eref: (float) AP1roG reference energy (default ``float('nan')``)

    :ecore: (float) core energy (default ``float('nan')``)

    :threshold: (float) optimization threshold for amplitudes (default ``1e-6``)

    :maxiter: (int) maximum number of iterations (default ``200``)

    :guess: (1-dim np.array) initial guess (default ``None``). In not provided, an initial guess containing random numbers in the interval :math:`(0,0.01]` is generated. The random guess preserves the symmetry of the PTb amplitudes, that is, :math:`t_{ij}^{ab}=t_{ji}^{ba}`. If a user-defined guess is provided, the elements of :math:`t_{ij}^{ab}` have to be indexed in C-like order

Example Python scripts
======================

MP2 calculation on the water molecule
-------------------------------------

This is a basic example on how to perform a RMP2 calculation in Horton. This script performs a RMP2 calculation on the water molecule using the cc-pVDZ basis set.

.. literalinclude:: ../data/examples/perturbation_theory/water_mp2.py
    :caption: data/examples/perturbation_theory/water_mp2.py
    :lines: 2-


PTa calculation on the water molecule
-------------------------------------

This is a basic example on how to perform a PTa calculation in Horton. This script performs a PTa calculation on the water molecule using the cc-pVDZ basis set.

.. literalinclude:: ../data/examples/perturbation_theory/water_pta.py
    :caption: data/examples/perturbation_theory/water_pta.py
    :lines: 2-


PTb calculation on the water molecule
-------------------------------------

This is a basic example on how to perform a PTb calculation in Horton. This script performs a PTb calculation on the water molecule using the cc-pVDZ basis set.

.. literalinclude:: ../data/examples/perturbation_theory/water_ptb.py
    :caption: data/examples/perturbation_theory/water_ptb.py
    :lines: 2-
