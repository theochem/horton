..
    : HORTON: Helpful Open-source Research TOol for N-fermion systems.
    : Copyright (C) 2011-2016 The HORTON Development Team
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

Release history
###############

Documentation of all HORTON versions can be found here: http://theochem.github.com/horton/

**June 17, 2016. Version 2.0.1**

   - QA framework for automatically testing pull requests on Github with Travis-CI.
   - Update config files for recent OSX, Ubuntu Linux and Fedora Linux versions.
   - Fix: contractions of Gaussian functions are normalized when creating new basis sets.
     (When loading wavefunctions from files that also contain a description of the basis
     set, the contractions are not renormalized for the sake of consistency.)
   - Fix: The numerical Poisson solver now also computes the correct asymptotics in the
     limit of small radii.
   - Fix: The two methods to project orbitals onto new basis sets (in
     ``horton/meanfield/project.py``) contained mistakes, which are now fixed.
   - Fix: Several unit tests using random data ocassionally failed, which is now fixed.
   - Many small bug fixes and corrections.


**June 11, 2015. Version 2.0.0**

   - The (orbital-optimized) AP1roG method (geminal-based wafecuntions).
   - Perturbation theory methods: MP2 (post-HF) and PTa and PTb (post-AP1roG).
   - Cholesky decomposition of the four-center integrals.
   - Installation instructions for Mac OS/X.
   - New file format: FCIDUMP.
   - Improvements in other file formats, e.g. all different conventions for the
     Molden file format are recognized automatically.
   - Model Hamiltonians.
   - Orbital localization (Pipek-Mezey).
   - Orbital entanglement analysis.
   - A Numerical Poisson solver (AIM analysis and pure KS-DFT implementation).
   - Evaluation of the kinetic energy density on a grid.
   - Projection of orbitals onto new basis sets.
   - Update to LibXC 2.2.2.
   - A lot of documentation, revamped website and code examples.
   - Many small cleanups and improvements under the hood.


**March 17, 2014. Version 1.2.1**

   - Update to LibXC-2.0.3 plus overview of the supported functionals in the
     documentation.
   - Update to h5py-2.2.1
   - Several bug fixes in the ESP fitting scripts.
   - Hu-Lu-Yang ESP cost function
   - Documentation for the ESP fitting scripts.
   - Mandatory output argument for most ``horton-*.py`` scripts.
   - Properly load fchk files from Gaussian calculations with Ghost atoms.
   - New script: ``horton-convert.py``. (Conversion between different file formats supported in HORTON.)
   - New script: ``horton-cubehead.py``. (Part of the ESP fitting scripts. A tool te generate economic grid specs for cubegen.)
   - Usability improvements in ``horton-atomdb.py``.
   - Skip expensive AIM computatoins by default in ``horton-wpart.py``
   - Documentation generation for C++ code with Doxygen and Breathe.
   - More covalent and van der Waals radii.
   - Several fixes in the CIF reader.
   - Improved EDIIS
   - Constructing a DFT/HF hamiltonian without Exchange term raises an error (unless idiot_proof is disable)
   - Additional basis sets
   - Several minor fixes and cleanups


**August 25, 2013. Version 1.2.0**

   - Gaussian/GAMESS wfn file reader. WFN files are now supported in
     ``horton-wpart.sh``. (Thanks to Farnaz!)
   - HORTON wavefunctions can now be written to the molden file format.
   - The efficiency of ``horton-wpart.sh`` has improved.
   - Added ``--lmax`` option to ``horton-wpart.sh`` and ``horton-cpart.sh`` to
     control the maximum angular momentum for the multipole analysis.
   - Fixed a division-by-zero-bug and a caching bug in the Iterative
     Stockholder scheme.
   - DIIS algorithms in ``horton.meanfield`` package: CDIIS [pulay1980]_, EDIIS
     and EDIIS+DIIS [kudin2002]_.
   - Improved efficiency of numerical integration in DFT hamiltonians.
   - A robust quadratic programming solver with linear (in)equality constraints.
     (This is used by EDIIS and Hirshfeld-E.)
   - Fix for compilation of libxc-2.0.2 with gfortran 4.8.1 and newer.
   - More detailed timer output. (Simplified usage of timer in source code.)
   - Improved screen output.
   - More documentation of the source code.
   - Several mistakes were fixed in the Gaussian basis set tutorial.
   - LineGrid and RectangleGrid for visualization purposes.
   - Various cleanups.


**July 22, 2013. Version 1.1.0**

   - Iterative Stockholder partitioning [lillestolen2008]_.
   - Pure (harmonic) multipoles in the AIM analysis.
   - Spin charges in the AIM analysis.
   - Switch to libxc-2.0.2.
   - New pruned atomic integration grids for elements H-La, Hf-Rn, with more
     levels of accuracy.
   - New radial integration grids with improved accuracy.
   - ADF is no longer supported in ``horton-atomdb.py``.
   - More efficient Becke weights.
   - Screen output and timer improvements.
   - A fast (approximate) evaluation of the electron density in
     ``horton-wpart.py``.
   - Many cleanups.


**July 5, 2013. Version 1.0.2**

   - Also support dynamic linking of libint and libx.
   - Switch to libint-2.0.3-stable.
   - Various cleanups.


**July 1, 2013. Version 1.0.1**

   - Two bug fixes related to reading Gaussian formatted checkpoint files.

     1. The Gaussian 03 FCHK format contains a spelling error ('independant'
        instead of 'independent'). This is fixed in Gaussian 09. Both variants
        are now properly handled by HORTON.
     2. Post-HF density matrices are read in properly.

   - Reorganization of mean-field code. It is now located in a sub package
     ``horton.meanfield``.
   - It is now impossible to start the SCF-ODA algorithm with a density matrix
     whose occupation numbers fall out of the admissible range. This prevents
     `fake` convergence to nonphysical solutions.
   - ESP fitting for isolated systems.


**May 23, 2013. Version 1.0**

   - This release mainly focuses on real-space density partitioning
     (atoms-in-molecules) methods.
   - Other major features include: import wavefunctions from various
     file formats, basic Hartree-Fock and DFT algorithms (making user of libint
     and libxc), pruned integration grids up to Ar, checkpointing, ...
   - Experimental features: ESP fitting of charges and related algorithms,
     currently only for 3D periodic systems.
