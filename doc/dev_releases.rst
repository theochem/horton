Release history
###############

* **June 22, 2013. Version 1.1.0**

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


* **June 5, 2013. Version 1.0.2**

   - Also support dynamic linking of libint and libx.
   - Switch to libint-2.0.3-stable.
   - Various cleanups.


* **June 1, 2013. Version 1.0.1**

   - Two bug fixes related to reading Gaussian formatted checkpoint files.

     1. The Gaussian 03 FCHK format contains a spelling error ('independant'
        instead of 'independent'). This is fixed in Gaussian 09. Both variants
        are now properly handled by Horton.
     2. Post-HF density matrices are read in properly.

   - Reorganization of mean-field code. It is now located in a sub package
     ``horton.meanfield``.
   - It is now impossible to start the SCF-ODA algorithm with a density matrix
     whose occupation numbers fall out of the admissible range. This prevents
     `fake` convergence to nonphysical solutions.
   - ESP fitting for isolated systems.


* **May 23, 2013. Version 1.0**

   - This release mainly focuses on real-space density partitioning
     (atoms-in-molecules) methods.
   - Other major features include: import wavefunctions from various
     file formats, basic Hartree-Fock and DFT algorithms (making user of libint
     and libxc), pruned integration grids up to Ar, checkpointing, ...
   - Experimental features: ESP fitting of charges and related algorithms,
     currently only for 3D periodic systems.
