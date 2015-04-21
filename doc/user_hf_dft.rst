.. _howtoscf:

How to use Horton as an Hartree-Fock/DFT program
################################################

Basic example
=============

This is a basic example of a Hartree-Fock computation in Horton. The input file
is just a small Python main program that uses the Horton library. The script
``data/examples/001_hf_water/run.py`` performs a HF/3-21G computation on water
and partitions the density with the Becke scheme:

.. literalinclude:: ../data/examples/001_hf_water/run.py

The molecular geometry is loaded from an XYZ file. It is also possible (yet less
convenient) to directly type the atomic coordinates in the python script.

In addition to this Python interface of Horton, several Python scripts were also
developed (``horton-*.py``) that are easier to use but only have a limited
functionality compared to the Python interface

TODO: more examples needed. Becke stuff should be removed.

How to add a custom external potential to a mean-field Hamiltonian
==================================================================

TODO

How to project orbitals on a new basis
======================================

TODO
