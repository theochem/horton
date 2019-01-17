..
    : HORTON: Helpful Open-source Research TOol for N-fermion systems.
    : Copyright (C) 2011-2019 The HORTON Development Team
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

How to use HORTON to convert between file formats?
##################################################

As discussed in the :ref:`using_horton_as_a_script` section, the (non-)optional
command-line arguments of the conversion scripts can be looked at by the ``--help`` flag:

.. code-block:: bash

    horton-convert.py --help
    horton-hdf2csv.py --help

``horton-convert.py`` -- Conversion between file formats supported by HORTON
============================================================================

The ``horton-convert.py`` script converts between various file formats supported by HORTON. To use it, simply provide the input and output files, as follows:

.. code-block:: bash

    horton-convert.py input_file output_file

The formats you are converting from/to are deduced from the input/output file name specified. This conversion only works if both file formats are recognized by HORTON, and the required information to write the output file is available in the input file. For more information, please refer to :ref:`ref_file_formats`.

This script is typically useful for archiving large amounts of data produced by other codes in text format. The following example converts a ``*.cube`` file to HORTON's internal HDF5 format, and then uses ``h5repack``, one of the HDF5 command-line tools, to compress the binary file.

.. code-block:: bash

    horton-convert.py density.cube density_cube.h5
    h5repack -v -f GZIP=1 density_cube.h5 density_cube_gzip.h5


.. _hdf2csv:

``horton-hdf2csv.py`` -- Conversion of HDF5 files to CSV format
===============================================================

The ``horton-hdf2csv.py`` script is used by providing the HDF5 and the CSV file names, respectively, as follows:

.. code-block:: bash

    horton-hdf2csv.py somefile.h5:path/in/hdf5/file otherfile.csv

For example, if ``horton-wpart.py`` was used to run an Extended Hirshfeld partitioning of a wavefunction in ``gaussian.fchk`` file format, the command below would convert the results into CSV file format:

.. code-block:: bash

    horton-hdf2csv.py wpart_he_results.h5 extended_hirshfeld.csv

This script was added for the convenience of those who are not familiar with the HDF5 files. If you are familiar with the HDF5 files, we strongly suggest you to directly process these files with conventional (Python) scripts. This would be far easier than interfacing with the CSV files. The `h5py <http://www.h5py.org/>`_ library is a great tool to make such custom scripts.
