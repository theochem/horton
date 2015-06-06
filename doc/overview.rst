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

To be completed soon ...


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
