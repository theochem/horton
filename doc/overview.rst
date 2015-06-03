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
  [marques2012]_ and `LibInt2 <https://github.com/evaleev/libint>`_.

* Support for many :ref:`ref_file_formats` to exchange data with other codes.

* The modular structure allows HORTON to be combined with custom developments.

**Electronic structure methods**

* Hamiltonians

    * Molecular electronic Hamiltonians in localized Gaussian basis sets

        * Four-center integrals computed with
          `LibInt2 <https://github.com/evaleev/libint>`_

        * :ref:`ref_gaussian_basis_standard_sets`

    * Model Hamiltonians

    * User-provided Hamiltonians

* Hatree-Fock and DFT methods:

    * Restricted and unrestricted orbitals

    * LDA, GGA and Hybrid GGA :ref:`ref_functionals`

    * Various SCF algorithms

* AP1roG

* Perturbation theory

**Post-processing**

* Several variants of the (iterative) Hirshfeld method

* Electrostatic Potential Fitting

* Orbital entanglement analysis

* Orbital localization
