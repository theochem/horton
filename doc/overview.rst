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

.. _overview:

HORTON Overview
###############

Our Manifesto
=============

HORTON is a *H*\ elpful *O*\ pen-source *R*\ esearch *TO*\ ol for *N*-fermion
systems, written primarily in the `Python language <https://www.python.org/>`_.
(HORTON is named for the `helpful pachyderm
<http://en.wikipedia.org/wiki/Horton_the_Elephant>`_, not the `Canadian caffeine
supply store <http://www.timhortons.com/>`_.) The ultimate goal of HORTON is to
provide a useful platform for testing new ideas related to the quantum many-body
problem, at reasonable computational cost. Although HORTON is primarily intended
as a `quantum-chemistry program
<http://en.wikipedia.org/wiki/List_of_quantum_chemistry_and_solid-state_physics_software>`_,
it can be used for model Hamiltonians and could be extended to nuclear physics.

What HORTON is (or hopes to be):

* HORTON is designed as a *helpful* framework for rapidly prototyping methods
  and testing ideas, together with utilities that help ensure that the resulting
  implementation is not too inefficient. HORTON is not designed to achieve
  bleeding-edge performance: readability, extensibility, and modifiability are
  often preferred over computational efficiency and code compactness when
  tradeoffs are essential.

* HORTON is, and always will be, *open source*, distributed under the :ref:`GNU
  General Public License <license_information>`. Its developers welcome new
  contributors.

* HORTON is a *research tool* for its contributors and users. Its functionality
  is naturally biased towards the interests of its developers, so its current
  focus is on low-cost *ab initio* electronic structure theory methods and
  post-processing tools for interpreting electronic structure calculations.
  Additional functionality will be provided by other developers and through
  interfaces to other programs. People interested in joining the HORTON
  development team should contact us through the `the HORTON mailing list
  <https://groups.google.com/forum/#!forum/horton-discuss>`_.

* HORTON can be used as a stand-alone program for the N-fermion problem using
  either input files or Python scripts. Although more work remains to be done,
  we hope that HORTON is also (1) easy to use as a procedure providing
  specialized functionality to other programs and (2) easy to use as a scripting
  language for processing results from other Schrödinger solvers (which may have
  alternative capabilities or better computational performance). To encourage
  others to use HORTON within their programs, we guarantee no major `API
  <http://en.wikipedia.org/wiki/Application_programming_interface>`_ changes for
  at least twelve months after each major release. Ensuring that HORTON can be
  used in all three ways (inside a larger program; as a stand-alone program;
  outside other programs) is a fundamental design principle of HORTON.


Main Features
=============

**Interoperability**

* Integration of the `LibXC
  <http://www.tddft.org/programs/octopus/wiki/index.php/Libxc>`_
  [marques2012]_ and `LibInt2 <https://github.com/evaleev/libint>`_
  [valeev2014]_.

* Support for many :ref:`ref_file_formats` to exchange data with other codes.

* The modular structure allows HORTON to be combined with custom developments.

**Electronic structure methods**

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

**Post-processing**

* Several variants of the (iterative) Hirshfeld method

* Electrostatic Potential Fitting

* Orbital entanglement analysis

* Orbital localization
