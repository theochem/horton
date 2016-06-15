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

.. _ref_grids:

Atomic integration grids
########################

Built-in grids
--------------

The table below lists the number of grid points in each built-in atomic
integration grid in HORTON. For the heavier elements (after Kr), the grids are
calibrated based on reference computations with pseudo potentials. The last
digit in the grid name refers to the number of significant digits obtained with
these grids on diatomic benchmarks.

.. include:: grids.rst.inc


.. _ref_grid_option:

The ``--grid`` option
---------------------

The scripts ``horton-atomdb.py`` and ``horton-wpart.py`` have a ``--grid``
option to specify the atomic integration grids. This option may be used to
control the numerical accuracy of the output. This option takes one argument
that may have any of the following four formats:

1. A filename of with atomic grid specifications in the same format as those
   found in ``${HORTONDATA}/grids``.

2. A filename of a file in ``${HORTONDATA}/grids`` without the extension
   ``.txt``. (These are discussed above.)

3. Equivalently to the previous point, one of the following names: ``coarse``,
   ``medium``, ``fine``, ``veryfine``, ``ultrafine`` and ``insane``.

4. A string with the following format: ``rname:rmin:rmax:nrad:nll``, with
   the following meaning for the keywords:

     * ``rname`` specifies the type of radial grid. It can be ``linear``,
       ``exp`` or ``power``.

     * ``rmin`` and ``rmax`` specify the first and the last radial grid point in
       angstroms.

     * ``nrad`` is the number of radial grid points.

     * ``nll`` is the number of points for the angular Lebedev-Laikov grid. It
       must be one of the following values: 6, 14, 26, 38, 50, 74, 86, 110, 146,
       170, 194, 230, 266, 302, 350, 434, 590, 770, 974, 1202, 1454, 1730, 2030,
       2354, 2702, 3074, 3470, 3890, 4334, 4802, 5294, 5810
