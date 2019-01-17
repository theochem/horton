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

.. _hamiltonian_io:

Dumping/Loading a Hamiltonian to/from a file
############################################

HORTON supports two formats for Hamiltonians: (i) an internal binary format based
on HDF5 (extension ``.h5``) and (ii) Molpro's FCIDUMP text format (containing ``FCIDUMP``
somewhere in the file name). The internal format is more flexible and can store
Hamiltonians in various ways. The FCIDUMP format is more restricted
but can be used to interface HORTON with different codes, e.g. Molpro. HORTON can
also load integrals from a Gaussian log file but this is **absolutely not**
recommended for any serious calculation.

For general information on how to load and dump data with HORTON in different
data file formats, refer to :ref:`ref_file_formats`.


.. _hamiltonian_io_internal:

HORTON's internal format
========================

Dumping
-------

You can store all of the operators in the atomic-orbital (AO) basis. For example,
you can store the core energy and all the integrals in the Cholesky decomposed format:

.. literalinclude :: ../data/examples/hamiltonian/dump_internal_ao.py
    :caption: data/examples/hamiltonian/dump_internal_ao.py
    :lines: 3-23

The internal format will store all attributes of the
:py:class:`~horton.io.iodata.IOData` instance. Note that the attributes
``coordinates``, ``numbers`` and ``title`` were loaded from the ``.xyz`` file
and will also be dumped into the internal format. Deleting an attribute, e.g.
``del mol.title``, or adding a new attribute, e.g. ``mol.egg = 'spam'``, will
result in the exclusion or the inclusion of the appropriate attribute,
respectively.

In the HDF5 file, all data is stored in binary form with full precision, which
means that all significant digits are written to file. In the example above, the
dumped HDF5 file will have the following layout:

.. code-block:: text

    $ h5dump -n hamiltonian_ao.h5
    HDF5 "hamiltonian.h5" {
    FILE_CONTENTS {
     group      /
     dataset    /coordinates
     dataset    /core_energy
     group      /er
     dataset    /er/array
     group      /kin
     dataset    /kin/array
     group      /lf
     group      /na
     dataset    /na/array
     dataset    /numbers
     dataset    /title
     }
    }

The attributes of the FCIDUMP format, ``one_mo``, ``two_mo``, ``core_energy``,
``nelec`` and ``ms2`` can also be used in the internal format. We can create an
empty :py:class:`~horton.io.iodata.IOData` instance and assign each of these
attributes. For example,

.. literalinclude::  ../data/examples/hamiltonian/dump_internal_ao_fcidump.py
    :caption: data/examples/hamiltonian/dump_internal_ao_fcidump.py
    :lines: 3-32

will results in the following HDF5 layout:

.. code-block:: text

    $ h5dump -n hamiltonian_ao_fcidump.h5
    HDF5 "hamiltonian_ao_fcidump.h5" {
    FILE_CONTENTS {
     group      /
     dataset    /core_energy
     group      /lf
     dataset    /ms2
     dataset    /nelec
     group      /one_mo
     dataset    /one_mo/array
     group      /two_mo
     dataset    /two_mo/array
     }
    }

Note that the integrals in this example are actually stored in the AO basis
(despite the ``_mo`` suffix). Read the section :ref:`user_hf_dft_preparing_posthf`
if you want to compute (and store) integrals in the molecular orbital (MO) basis.

Loading
-------

You can load integrals from the HORTON's internal format, as follows:

.. literalinclude :: ../data/examples/hamiltonian/load_internal_ao.py
    :caption: data/examples/hamiltonian/load_internal_ao.py
    :lines: 3-11


.. _hamiltonian_io_fcidump:


FCIDUMP format
==============

Dumping
-------

The FCIDUMP format is useful when exchanging Hamiltonians with different codes.
Unlike the internal format, there are some restrictions:

1. One-body operators (usually kinetic energy and nuclear attraction)
   must all be added into a single two-index object.
2. Integrals can only be stored in a restricted (MO) basis set, i.e. using
   different basis sets for the alpha and beta orbitals is not possible.
3. Two-electron integrals must be stored as a ``DenseFourIndex`` object, so
   the Cholesky decomposition of the ERI is not supported.

The FCIDUMP format is normally used for storing integrals in the MO basis
but in the following example we will use the AO basis. Read the section
:ref:`user_hf_dft_preparing_posthf` if you want to compute (and store) integrals
in the MO basis. The storage of integrals in AO basis in FCIDUMP format is as follows:

.. literalinclude :: ../data/examples/hamiltonian/dump_fcidump_ao.py
    :caption: data/examples/hamiltonian/dump_fcidump_ao.py
    :lines: 3-28

In this example, we set the :py:class:`~horton.io.iodata.IOData` attributes by
using keyword arguments in the constructor. The file ``hamiltonian_ao.FCIDUMP``
will contain the following:

.. code-block:: text

     &FCI NORB=28,NELEC=20,MS2=0,
      ORBSYM= 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
      ISYM=1
     &END
     5.9786161694613265e+00    1    1    1    1
    -1.0433683639500952e+00    2    1    1    1
     1.9823368119430076e-01    2    1    2    1
     6.1978840317716133e-01    2    2    1    1
    .
    .
    .
    -1.4611205008249136e+01   27   27    0    0
    -8.5804159203863803e-12   28   14    0    0
    -1.4611205008249140e+01   28   28    0    0
     2.0000000000000000e+01    0    0    0    0

This file is divided into two blocks. The first block (between ``&FCI`` and
``&END``) contains system-specific information:

    :NORB: number of orbitals/basis functions
    :NELEC: number of electrons
    :MS2: spin projection
    :ORBSYM: irreducible representation of each orbital
    :ISYM: total symmetry of the wavefunction

The second block (after ``&END``) contains the one- and two-electron integrals
and the core energy in that order:

1. All symmetry-unique elements of the two-electron integrals are listed,
   where the first column is the value of the integral, followed by the orbital
   indices. Note that the orbital indices start from one and that the orbital indices
   (``i j k l``) in an FCIDUMP file are written in chemists' notation,

  .. math::

      (ij\vert kl) = \langle ik \vert jl \rangle = \int \phi_i^*(\mathbf{x}_1)
      \phi_k^*(\mathbf{x}_2) \frac{1}{r_{12}} \phi_j(\mathbf{x}_1) \phi_l(\mathbf{x}_2)
      d\mathbf{x}_1 d\mathbf{x}_2

2. All symmetry-unique elements of the one-electron integrals are listed, where
   the first column is the value of the integral, followed by the orbital
   indices. Note that the last two columns contain zeros.

3. Core energy (for instance, the nuclear repulsion term, etc.) is written
   on the last line with all orbital indices equal 0.

If the value of an integral is zero, the corresponding line is not included in
the the FCIDUMP file. It's important to note that HORTON does not (yet) support
geometric symmetries, so all the orbitals will be assigned a symmetry label 1,
which corresponds to the point group C1 in the FCIDUMP format.

Loading
-------

You can load integrals, stored in an FCIDUMP file, as follows:

.. literalinclude :: ../data/examples/hamiltonian/load_fcidump_ao.py
    :caption: data/examples/hamiltonian/load_fcidump_ao.py
    :lines: 3-12
