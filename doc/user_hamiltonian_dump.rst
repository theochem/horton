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

.. _hamiltonian_dump:

Dumping a Hamiltonian to a file
###############################

Horton supports two formats for storing a Hamiltonian in a file: (i) an internal binary format based on HDF5 (extension ``.h5``) and Molpro's FCIDUMP text format (containing ``FCIDUMP`` somewhere in the file name). The internal format is more flexible and can store a Hamiltonian on various different ways. The FCIDUMP format as more restrictive but can be used to interoperate with different codes, e.g. Molpro.

All input and output of data in Horton is done through the :py:class:`horton.io.iodata.IOData` class. Dumping data to a file takes the following three steps:

1. Create an instance of the ``IOData`` class (from scratch, or by loading it from a file),
2. Sets some attributes,
3. Call the :py:meth:`~horton.io.iodata.IOData.to_file` method to actually dump the info to the file.

One can also combine steps 1 and 2 by passing the attributes of interest to the constructor. These steps will be illustrated in the examples below.

The list of attributes recognized by (some) output formats is documented here: :py:class:`horton.io.iodata.IOData`. The method :py:meth:`~horton.io.iodata.IOData.to_file` just takes the filename as argument. The filename is used to determine the file format. When it has the extension ``.h5``, the internal format is used. When the file contains ``FCIDUMP``, the FCIDUMP format is used.

.. _hamiltonian_dump_internal:

Horton's Internal format
========================

One can store all the separate operators in the atomic-orbital (AO) basis, here using a Cholesky decomposition of the four-center integrals, as follows:

.. literalinclude :: ../data/examples/hamiltonian/dump_internal_ao.py
    :caption: data/examples/hamiltonian/dump_internal_ao.py
    :lines: 2-

Note that the attributes ``coordinates``, ``numbers`` and ``title``, which were loaded from the ``.xyz`` file, will also be dumped in the internal format, unless you explicitly remove them first with e.g. ``del mol.title``. The internal format will just store any attribute of the ``IOData`` object, not just the ones that are documented, see :py:class:`horton.io.iodata.IOData`. So, you may assign other attributes as well, e.g. ``obasis`` or ``olp``, if that is convenient.

In the HDF5 file, all data is stored binary in full precision. The layout of the
HDF5 file in this example is as follows:

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
     group      /na
     dataset    /na/array
     dataset    /numbers
     dataset    /title
     }
    }

The attributes used for the FCIDUMP format, ``one_mo``, ``two_mo``, ``core_energy``, ``nelec`` and ``ms2`` can also be used for the internal format. The following example shows how this can be done with a ``IOData`` object created from scratch.

.. literalinclude::  ../data/examples/hamiltonian/dump_internal_ao_fcidump.py
    :caption: data/examples/hamiltonian/dump_internal_ao_fcidump.py
    :lines: 2-

which results in the following HD5 layout:

.. code-block:: text

    $ h5dump -n hamiltonian_ao_fcidump.h5
    HDF5 "hamiltonian_ao_fcidump.h5" {
    FILE_CONTENTS {
     group      /
     dataset    /core_energy
     dataset    /ms2
     dataset    /nelec
     group      /one_mo
     dataset    /one_mo/array
     group      /two_mo
     dataset    /two_mo/array
     }
    }

Note that the integrals in this example are actually stored in the AO basis (unlike the ``_mo`` suffix suggests). Read the section :ref:`user_hf_dft_preparing_posthf` if you want to compute (and store) integrals in the molecular-orbital (MO) basis.


.. _hamiltonian_dump_fcidump:

FCIDUMP format
==============

The FCIDUMP format is mainly useful for exchanging Hamiltonians with different codes. Compared to the internal format, there are some restrictions:

1. The one-body terms must all be added into a single operator.
2. The integrals can only be stored in a restricted (MO) basis set, i.e. so no different basis sets for the alpha and beta orbitals are possible.
3. The two-electron integrals must be stored in a ``DenseFourIndex`` object, so the Cholesky decomposition of the ERI is not supported.

The FCIDUMP format is normally only used for storing integrals in the MO basis but the example here will only consider the AO basis. Read the section :ref:`user_hf_dft_preparing_posthf` if you want to compute (and store) integrals in the molecular-orbital (MO) basis. The usage is as follows:

.. literalinclude :: ../data/examples/hamiltonian/dump_fcidump_ao.py
    :caption: data/examples/hamiltonian/dump_fcidump_ao.py
    :lines: 2-

This example shows how the ``IOData`` attributes can be set by giving keyword arguments to the constructor. The file ``hamiltonian_ao.FCIDUMP`` will contain the following:

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

The input file is divided into two blocks. The first block (between ``&FCI`` and ``&END``) contains system-specific information:

    :NORB: number of orbitals/basis functions
    :NELEC: number of electrons
    :MS2: spin projection
    :ORBSYM: irreducible representation of each orbital
    :ISYM: total symmetry of the wavefunction

The second block (after ``&END``) contains the one- and two-electron integrals as well as the core energy:

* First, all symmetry-unique elements of the two-electron integrals are listed, where the first column is the value of the integral, followed by the orbital indices. (Orbital indices start counting from one.) Note that the orbital indices (``i j k l``) in an FCIDUMP file are written in chemists' notation,

  .. math::

      (ij\vert kl) = \langle ik \vert jl \rangle = \int \phi_i^*(\mathbf{x}_1) \phi_k^*(\mathbf{x}_2) \frac{1}{r_{12}} \phi_j(\mathbf{x}_1) \phi_l(\mathbf{x}_2) d\mathbf{x}_1 d\mathbf{x}_2

* Second, all symmetry-unique elements of the one-electron integrals are listed, where again the first column is the value of the integral, followed by the orbital indices. Note that the last two columns contain zeros.

* Finally, the core energy (for instance, the nuclear repulsion term, etc.) is written on the last line with all orbital indices equal 0.

If any of the value of an integral is zero, the corresponding line is not included in the the FCIDUMP file.
