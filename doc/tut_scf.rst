How to use Horton as an SCF program
###################################

Note that this part of Horton is very experimental and will undergo a
serious redesign in the near future.

Basic example
=============

This is a basic example of a Hartree-Fock computation in Horton. The input file
is just a small Python main program that uses the Horton library. The script
``examples/001_hf_water/run.py`` performs a HF/3-21G computation on water and
partitions the density with the Becke scheme:

.. literalinclude:: ../examples/001_hf_water/run.py

The molecular geometry is loaded from an XYZ file. It is also possible (yet less
convenient) to directly type the atomic coordinates in the python script.

In addition to this Python interface of Horton, several Python scripts were also
developed (``horton-*.py``) that are easier to use but only have a limited
functionality compared to the Python interface
