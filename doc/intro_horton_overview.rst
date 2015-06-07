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

.. _intro_horton_overview:

HORTON Overview
###############

Our Manifesto
=============

HORTON is a *H*\ elpful *O*\ pen-source *R*\ esearch *TO*\ ol for *N*-fermion
systems, written primarily in the `Python programming language <https://www.python.org/>`_.
(HORTON is named after the `helpful pachyderm
<http://en.wikipedia.org/wiki/Horton_the_Elephant>`_, not the `Canadian caffeine
supply store <http://www.timhortons.com/>`_.) The ultimate goal of HORTON is to
provide a platform for testing new ideas on the quantum many-body
problem at a reasonable computational cost. Although HORTON is primarily designed
to be a `quantum-chemistry program
<http://en.wikipedia.org/wiki/List_of_quantum_chemistry_and_solid-state_physics_software>`_,
it can perform computations involving model Hamiltonians, and could be extended for computations in nuclear physics.

What HORTON is (or is hoped to be):

* HORTON is designed to be a *helpful* framework for rapidly prototyping methods
  and testing ideas, together with utilities that ensure the resulting
  implementation is not too inefficient. HORTON is not designed to achieve
  bleeding-edge performance: readability, extensibility, and modifiability are
  often preferred over computational efficiency and code compactness when
  tradeoffs are essential.

* HORTON is, and always will be, *open source*, distributed under the :ref:`GNU
  General Public License <license_information>`. The HORTON development team always welcomes new
  contributions.

* HORTON is a *research tool* for both its developers and users. As a result, the available functionality
  is naturally biased towards the interests of its developers. The current
  focus of HORTON is on low-cost *ab initio* electronic structure theory methods and
  post-processing tools for interpreting electronic structure calculations.
  Additional functionality can be provided by other developers, and through
  interfaces to various programs. If you are interested in joining the HORTON
  development team, please contact us through the `the HORTON mailing list
  <https://groups.google.com/forum/#!forum/horton-discuss>`_.

* HORTON can be used as a stand-alone program for the *N-fermion* problems using
  either input files or Python scripts. Although more work remains to be done,
  we hope that HORTON is easy to use as (1) a procedure providing
  specialized functionality to other programs and (2) a scripting
  language for processing results from other Schrödinger solver software (which may have
  alternative capabilities or better computational performance). To faciliate
  the use of HORTON within other programs, we guarantee no major `API
  <http://en.wikipedia.org/wiki/Application_programming_interface>`_ changes for
  at least twelve months after each major release. Ensuring that HORTON can be
  used in all three ways (as part of a larger program; as a stand-alone program;
  outside other programs) is a fundamental design principle of HORTON.


Main Features
=============

**Interoperability**

* Integrating the `LibXC
  <http://www.tddft.org/programs/octopus/wiki/index.php/Libxc>`_
  [marques2012]_ and `LibInt2 <https://github.com/evaleev/libint>`_
  [valeev2014]_ libraries.

* Supporting many :ref:`ref_file_formats` to exchange data with other codes.

* Adapting a modular structure to easily combine HORTON with custom developments.

**Electronic Structure Methods**

* Hamiltonians

    * Molecular electronic Hamiltonians in localized Gaussian basis sets

        * Four-center integrals computed with
          `LibInt2 <https://github.com/evaleev/libint>`_ [valeev2014]_

        * :ref:`ref_gaussian_basis_standard_sets`

    * Model Hamiltonians

    * User-provided Hamiltonians

* Hartree-Fock and DFT methods:

    * Restricted and unrestricted orbitals

    * LDA, GGA and Hybrid GGA :ref:`ref_functionals`

    * Various SCF algorithms

* Geminals-based methods. Specifically, Antisymmetry Product of 1
  reference-orbital Geminals (AP1roG)

* Perturbation theory

**Post-Processing**

* Several variants of the Hirshfeld atoms-in-molecules partitioning scheme

* Electrostatic potential fitting of atomic charges

* Orbital entanglement analysis

* Orbital localization
