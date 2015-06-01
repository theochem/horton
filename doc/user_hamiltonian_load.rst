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

.. _hamiltonian_load:

Loading a Hamiltonian from a file
#################################

Loading a Hamiltonian from a file is just the reverse process of :ref:`hamiltonian_dump`. For loading, the same two formats are supported: (i) an internal binary format based on HDF5 (extension ``.h5``) and Molpro's FCIDUMP text format (containing ``FCIDUMP`` somewhere in the file name).

For a general information on how to load and dump data with Horton in different data file formats, refer to :ref:`ref_file_formats`.


Horton's Internal format
========================

For more details in Horton's internal format, read :ref:`hamiltonian_dump_internal`.
The usage pattern is as follows:

.. literalinclude :: ../data/examples/hamiltonian/load_internal_ao.py
    :caption: data/examples/hamiltonian/load_internal_ao.py
    :lines: 2-


FCIDUMP format
==============

For more details in Horton's internal format, read :ref:`hamiltonian_dump_fcidump`.
The usage pattern is as follows:

.. literalinclude :: ../data/examples/hamiltonian/load_fcidump_ao.py
    :caption: data/examples/hamiltonian/load_fcidump_ao.py
    :lines: 2-
