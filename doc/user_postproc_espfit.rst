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

Electrostatic potential fitting
###############################


.. _user_espfit_introduction:

Introduction
============

The general idea behind the electrostatic potential (ESP) fitting tools in HORTON is to allow:

 * ESP-fitting: Fitting atomic charges to reproduce the actual molecular ESP, and

 * ESP-testing: Assessing how well a given set of charges reproduce the actual molecular ESP

in isolated and 3D-periodic systems with compatible
methodologies. These procedures are broken into small steps, such
that you can carry out more specialized ESP fittings.

These two computations, fitting charges and testing the quality of given charges
to reproduce the ESP, can be carried out once an ESP cost function
is constructed. In HORTON, this ESP cost function takes the following form:

.. math::
    \text{COST}(\{q\}, \Delta V_\text{ref}) = \int_V d\mathbf{r} w(\mathbf{r}) \left(V_\text{ai}(\mathbf{r}) - \sum_A q_A V_\text{point}(\mathbf{r} - \mathbf{R}_A) - \Delta V_\text{ref} \right)^2

where the symbols have the following meaning:

* :math:`\{q\}`: the set of the atomic charges
* :math:`V`: the volume where the point-charge ESP (approximated molecular ESP) is compared to the
  reference ESP (actual molecular ESP). In isolated systems, this is some volume surrounding the molecule.
  In periodic systems, this is the volume of a single primitive cell.
* :math:`w(\mathbf{r})`: the weight function with range [0, 1] that selects the
  relevant parts of the volume for the ESP fit. In HORTON, the weight function
  of Hu, Lu and Yang is used. [hu2007]_ This weight function is zero far away
  from and inside the atoms. It becomes one in the part of the electron density
  tail where non-bonding contacts are typically found. The weight function
  smoothly varies from 0 to 1. In the Hu-Lu-Yang paper, a pro-density is
  recommended for this purpose while we recommend the use of a proper all-electron
  density or a corrected pseudo-density.
* :math:`V_\text{ai}(\mathbf{r})`: the actual molecular ESP from `ab initio` or DFT calculations obtained
  from programs, like Gaussian, CP2K, VASP, etc.
* :math:`q_A`: the charge of atom :math:`A`.
* :math:`V_\text{point}(\mathbf{r} - \mathbf{R}_A)`: the ESP due to a unit point
  charge. In an isolated molecule, this is simply :math:`1/|\mathbf{r} - \mathbf{R}_A|`
  (in atomic units). For periodic systems, this part is quite a bit more involved.
* :math:`\Delta V_\text{ref}`: constant that may account for differences in
  reference between the ab initio ESP and the point charge ESP.
  The need for such a term was established in the REPEAT paper for ESP fitting
  in 3D periodic systems. [campana2009]_

This approach differs from traditional ESP fitting methods (RESP, MKS, CHELPG)
in the sense that the cost function is defined as an integral rather than a sum
over sample points. This idea comes from the Hu-Lu-Yang paper. [hu2007]_
The main difference with the REPEAT method [campana2009]_ is that the selected
volume for the fit is defined by a smooth weight function. When this weight
function is based on an all-electron density, there is no need to define atomic
radii as is done usually in ESP fitting.

The cost function is obviously a quadratic function of the unknown parameters,
:math:`X`, the charges (and :math:`\Delta V_\text{ref}` in the case of
a 3D periodic system), which we can always write as follows:

.. math::
    \text{COST}(X) = X^T A X - 2 B^T X - C

The script ``horton-esp-cost.py`` constructs the matrix :math:`A`, the vector
:math:`B` and the constant :math:`C`. These results are stored in an HDF5 file
that is subsequently used by the ``horton-esp-fit.py`` script to obtain
ESP-fitted charges, or by the ``horton-esp-test.py`` script to evaluate the
cost function for a given set of charges.

In the output of these three scripts, the following quantities are reported (if
applicable):

* :math:`\text{RMSD} = \sqrt{\text{COST}/V_w}`
* :math:`\text{RRMSD} = \sqrt{\text{COST}/C}` (expressed as a percentage)

where the cost is computed with a certain set of charges (and the optimal
:math:`\Delta V_\text{ref}` when applicable). The parameter :math:`V_w` is the
integral of the weight function. When these charges are the ESP-fitted charges,
:math:`\text{COST}`, :math:`\text{RMSD}` and :math:`\text{RRMSD}` are minimal.
When the charges are all set to zero, you obtain the worst-case value for
:math:`\text{COST}=C`, :math:`\text{RMSD}=\sqrt{C/V_w}` and
:math:`\text{RRMSD}=100\%`. (You can still obtain a higher
:math:`\text{COST}`, but only when using completely insensible charges. When that is
the case, it is better to set the charges equal to zero to model the ESP.)


``horton-esp-cost.py`` -- Set up an ESP cost function
=====================================================

The ``horton-esp-cost.py`` script can be used to construct the cost function for
the ESP fitting or charge testing.

For an all-electron calculation, it is recommended to use:

.. code-block:: bash

    horton-esp-cost.py esp.cube cost.h5 --wdens=rho.cube --pbc={000|111}

where ``esp.cube`` is a Gaussian cube file containing the actual molecular ESP on a grid. The results will
be written in ``cost.h5``. The option ``--wdens=rho.cube`` implies that the
weight function is constructed according to the Hu-Lu-Yang method [hu2007]_ using the given
density cube file. The option ``--pbc`` can be used to construct a cost
function for an isolated system (``000``) or a 3D periodic system (``111``).

When a pseudo-potential computation is used, the density cube file contains
regions of low electron density close to the nucleus. These regions may not be
excluded from the fit with the Hu-Lu-Yang weight function. Therefore, HORTON
allows you to build the weight function as a product of several factors:
:math:`w(\mathbf{r}) = w_1(\mathbf{r})w_2(\mathbf{r})w_3(\mathbf{r})w_4(\mathbf{r}) \ldots`, where the
first one is typically the Hu-Lu-Yang weight function and additional weight
functions can be included for every pseudo-core,

For a pseudo-potential calculation, it is recommended to use:

.. code-block:: bash

    horton-esp-cost.py esp.cube cost.h5 --wdens=rho.cube --pbc={000|111} \
                       --wnear Z1:r1:gamma1 [Z2:r2:gamma2 ...]

where the new option, ``--wnear``, takes at least one
argument. This argument must be presented for every element in the system for
which a pseudo potential is used. The parameter ``Z1`` is the atomic number of atom 1. The
parameter ``r1`` is the radius of the core region that you would like to
exclude for atom 1. The parameter ``gamma1`` determines how quickly the weight factor for
elements ``Z1`` switches from 0 (inside a sphere with radius ``r1``) to 1
(outside a sphere with radius ``r1``). Both ``r1`` and ``gamma1`` must be given
in angstrom.

The ``horton-esp-cost.py`` script has several more options. As discussed in the
:ref:`using_horton_as_a_script` section, the commands below describe
the arguments that each script take:

.. code-block:: bash

    horton-esp-cost.py --help
    horton-esp-fit.py  --help
    horton-esp-test.py --help
    horton-cubehead.py --help

``horton-esp-fit.py`` -- Fit atomic charges to the ESP
======================================================

Once a cost function is constructed, it can be used to estimate the atomic charges by
minimizing the cost function. A bare-bones fit can be carried out as follows:

.. code-block:: bash

    horton-esp-fit.py cost.h5 charges.h5

where the file ``cost.h5`` is the file generated by the
``horton-esp-cost.py`` script.

Useful ESP-fitted charges typically involve a more advanced minimizations of
the ESP cost function, for example, by adding constraints, restraints,
transforming to bond-charge increments, fitting to several different molecules
concurrently, etc. Such advanced features are not supported in
``horton-esp-fit.py`` script, but you can implement these in customized scripts
that use the the ESP cost functions obtained from the ``horton-esp-cost.py`` script.


``horton-esp-test.py`` -- Test the ESP quality of a given set of charges
========================================================================

The ``horton-esp-test.py`` script can be used to test the quality of a given set of
charges to reproduce the ESP. These charges are typically obtained
with ``horton-wpart.py``. This script can be used as follows:

.. code-block:: bash

    horton-esp-test.py cost.h5 wpart.h5:hi/charges wpart_espcost.h5:hi

The first file, ``cost.h5``, contains the cost function and is generated with the
``horton-esp-cost.py`` script. The second file, ``wpart.h5``, contains the given atomic chanrges
and, for example, is generated with ``horton-wpart.py gaussian.fchk wpart.h5:hi hi atoms.h5``. The last file,
``wpart_espcost.h5``, will contain the assessment results in the HDF5 group ``hi``.


Making nice cube files with Gaussian
====================================

HORTON contains an auxiliary tool, ``horton-cubehead.py``, to prepare an input
header for a cube file aligned with the molecule of interest. This is more
efficient than the default settings of cubegen, which makes a significant difference in
disk space when working with molecular databases. For occasional use,
``horton-cubehead.py`` is probably an overkill. The script is used as follows:

.. code-block:: bash

    horton-cubehead.py structure.xyz cubehead.txt

The ``cubehead.txt`` file will contain something along the following lines:

.. code-block:: text

    0   16.5695742234   -2.4411573645  -11.3378429796
  -61   -0.0000100512    0.0000288090    0.3779452256
   61   -0.2210334948    0.3065726468   -0.0000292468
   65   -0.3065726480   -0.2210334949    0.0000086952

This file can be used for the cubegen utility as follows:

.. code-block:: bash

    cubegen 0 fdensity=scf somefile.fchk rho.cube -1 < cubehead.txt
    cubegen 0 potential=scf somefile.fchk esp.cube -1 < cubehead.txt

where ``scf`` must be replaced by the type of wavefunction to be analyzed. Read
the `cubegen manual <http://www.gaussian.com/g_tech/g_ur/u_cubegen.htm>`_ for
more details.


Making cube files with CP2K
===========================

CP2K can generate cube files for periodic systems. You have to use the input sections `E_DENSITY_CUBE <http://manual.cp2k.org/trunk/CP2K_INPUT/FORCE_EVAL/DFT/PRINT/E_DENSITY_CUBE.html>`_ and `V_HARTREE_CUBE <http://manual.cp2k.org/trunk/CP2K_INPUT/FORCE_EVAL/DFT/PRINT/V_HARTREE_CUBE.html>`_ for density and ESP cube files, respectively. (The name ``V_HARTREE_CUBE`` is misleading. CP2K will print out the potential due to the electrons `and` the nuclei felt by a point charge with charge :math:`-e`.)


Python interface to the ESP fitting code
========================================

You can use the ESP cost function constructed by ``horton-esp-cost.py`` script to
implement customized charge fitting protocols, e.g. using bond-charge
increments, constraints or hyperbolic restraints. At the beginning of such
a customized script, the cost function can be loaded as follows:

.. code-block:: python

    cost = load_h5("cost.h5")['cost']

The ``cost`` object is an instance of the
:py:class:`~horton.espfit.cost.ESPCost` class. This
instance can, for example, be used to evaluate the ESP cost or its gradient for a
given array of atomic charges:

.. code-block:: python

    print cost.value(charges)
    print cost.gradient(charges)

If desired, you can also directly access :math:`A`, :math:`B`, :math:`C` that
define the quadratic cost functions:

.. code-block:: python

    print cost._A
    print cost._B
    print cost._C

For more information, please refer to :ref:`user_espfit_introduction`.
