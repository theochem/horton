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

.. _user_other_numerical_integration:

How to use HORTON for numerical integration?
############################################

.. _user_other_numerical_integration_specify:

Building an integration grid
============================

HORTON primarily makes use of Becke-Lebedev grids for numerical integrations over
molecular volumes. It can automatically construct Becke-Lebedev integration grids for a given
molecular geometry. This is needed for DFT computations (typically to
evaluate the exchange-correlation functional) or atoms-in-molecules
analyses. To familiarize yourself with these grids, please refer to [becke1988_multicenter]_ and
[lebedev1999]_.

In a nutshell, a Becke-Lebedev integration grid works as follows. Say, you are interested in
computing the integral of :math:`f(\mathbf{r})` over a molecular volume. In
practice, the integrand is often a functional of the electron density, so similar to the electron
density, it contains sharp spikes close to the atomic nuclei. In order to properly integrate all these
unsmooth spikes, the molecular integral is first split into atomic
contributions:

.. math::
    \int f(\mathbf{r}) d\mathbf{r} = \sum_A \int w_A(\mathbf{r}) f(\mathbf{r}) d\mathbf{r}

where :math:`w_A(\mathbf{r})` is the atomic weight function for atom A. This function is
1 close to the nucleus of atom A, and goes to zero upon leaving the domain of atom A. Every
atomic integral is then computed on a grid in a spherical coordinate system centered on
the nucleus of the corresponding atom. This atomic grid is typically a product grid, where different
one-dimensional radial grids are possible, and the Lebedev-Laikov grids are always used
for the angular part. Putting all these together, you can always approximate the numerical
integration as follows:

.. math::
    \int f(\mathbf{r}) d\mathbf{r} \approx \sum_{i=1}^{N_\text{grid}} w_i f(\mathbf{r}_i)

where :math:`N_\text{grid}` is the number of grid points, :math:`w_i` are the
integration grid weights and :math:`\mathbf{r}_i` are the integration grid
points.

The default Becke-Lebedev grid is constructed by specifying three non-optional arguments, as follows:

.. code-block:: python

    grid = BeckeMolGrid(coordinates, numbers, pseudo_numbers)

where ``coordinates`` is an array containing the Cartesian coordinates of the atoms,
``numbers`` is an array containing the (integer) atomic numbers of the atoms, and ``pseudo_numbers``
is an array containing the (float) effective core charges of the atoms. These arrays are always
available when you load a molecule from a file. For example:

.. code-block:: python

    mol = IOData.from_file('water.xyz')
    grid = BeckeMolGrid(mol.coordinates, mol.numbers, mol.pseudo_numbers)

The fourth optional argument can be provided to control the accuracy of the integration
grid, e.g.:

.. code-block:: python

    grid = BeckeMolGrid(coordinates, numbers, pseudo_numbers, 'fine')

The available levels of accuracy for the built-in numerical integration grids are
documented in :ref:`ref_grids`. You can also control in more detail the
radial and angular components of the integration grids. For more details,
please refer to the API documentation of :py:class:`horton.grid.molgrid.BeckeMolGrid`
and :py:class:`horton.grid.atgrid.AtomicGridSpec`.


Computing an integral involving the electron density
====================================================

This section assumes that the following objects are already available:

* ``obasis``: an orbital basis set object
* ``dm_full``: a spin-summed density matrix
* ``grid``: a Becke-Lebedev integration grid as introduced above.

If you are not familiar with the ``obasis`` and ``dm_full`` objects, please
refer to :ref:`user_molecularham_basis` and :ref:`user_hf_dft`,
respectively. Note that the density matrix can also be loaded from a file
instead of being computed by HORTON. This is demonstrated in an example at the
end of this section.

First, you must evaluate the integrand on all the points of the integration grid.
In case of the electron density, this can be done as follows:

.. code-block:: python

    rho = obasis.compute_grid_density_dm(dm_full, grid.points)

It is assumed that ``dm_full`` is symmetric. Several other quantities can also
be evaluated on the grid. For more details, please refer to:

* :py:meth:`~horton.gbasis.cext.GOBasis.compute_grid_density_dm`
* :py:meth:`~horton.gbasis.cext.GOBasis.compute_grid_gradient_dm`
* :py:meth:`~horton.gbasis.cext.GOBasis.compute_grid_gga_dm`
* :py:meth:`~horton.gbasis.cext.GOBasis.compute_grid_kinetic_dm`
* :py:meth:`~horton.gbasis.cext.GOBasis.compute_grid_hessian_dm`
* :py:meth:`~horton.gbasis.cext.GOBasis.compute_grid_mgga_dm`
* :py:meth:`~horton.gbasis.cext.GOBasis.compute_grid_hartree_dm`
* :py:meth:`~horton.gbasis.cext.GOBasis.compute_grid_esp_dm`
* :py:meth:`~horton.gbasis.cext.GOBasis.compute_grid_orbitals_exp`

Integrating the electron density by itself results in the total number of electrons.
This is a simple way to verify the accuracy of the integration grid.

.. code-block:: python

    print grid.integrate(rho)

Since ``rho`` is simply a Numpy array, it can be manipulated easily to compute
functionals of the electron density, e.g.

.. code-block:: python

    print grid.integrate(rho**(4.0/3.0))

You can also use the ``grid.points`` array to evaluate other expectation values
numerically, e.g. the following snippet evaluates the expectation value of
:math:`\vert\mathbf{r}\vert=(x^{2}+y^{2}+z^{2})^{0.5}`:

.. code-block:: python

    r = (grid.points[:,0]**2 + grid.points[:,1]**2 + grid.points[:,2]**2)**0.5
    grid.integrate(rho, r)

As shown in the above snippet, the ``integrate`` method can take multiple
one-dimensional arrays that are all multiplied before integration.

The following script is a complete example for computing the expectation value
of :math:`\vert\mathbf{r}\vert=(x^{2}+y^{2}+z^{2})^{0.5}`
for a molecular wave-function loaded from a file.

.. literalinclude:: ../data/examples/grid/expectation_r.py
    :caption: data/examples/grid/expectation_r.py
    :lines: 3-24


Constructing a one-body operator from a real-space potential
============================================================

This section assumes that the following objects are already available:

* ``obasis``: an orbital basis set object
* ``lf``: an instance of ``DenseLinalgFactory`` or ``CholeskyLinalgFactory``
* ``dm_full``: a spin-summed density matrix
* ``grid``: a Becke-Lebedev integration grid as introduced above.

If you are not familiar with the ``obasis`` or ``lf`` object, please refer to
:ref:`user_molecularham_basis`. The density matrix can either be read from a
file or computed with HORTON. For more information, please refer to :ref:`user_hf_dft`.

Given a multiplicative potential, its expectation value is written as:

.. math::

    \langle V \rangle = \int \rho(\mathbf{r}) V(\mathbf{r}) d\mathbf{r}.

Expanding the orbitals in a local basis set results in:

.. math::

    \langle V \rangle = \sum_{\mu\nu} D_{\mu\nu} \mathcal{V}_{\nu\mu}

where :math:`D_{\mu\nu}` is the spin-summed density matrix. The matrix
:math:`\mathcal{V}_{\nu\mu}` is defined as

.. math::

    \mathcal{V}_{\nu\mu} = \int V(\mathbf{r}) b_\nu^*(\mathbf{r}) b_\mu(\mathbf{r}) d\mathbf{r}

where :math:`b_\mu(\mathbf{r})` are the orbital basis functions. Such matrices
can be constructed with the
:py:meth:`~horton.gbasis.cext.GOBasis.compute_grid_density_fock` method. This method is
also useful when applying the chain rule to construct the contribution of a density functional
to a Fock matrix:

.. math::

    \frac{\partial E[\rho]}{\partial D_{\nu\mu}} = \int \frac{\delta E[\rho]}{\delta \rho(\mathbf{r})}  b_\nu^*(\mathbf{r}) b_\mu(\mathbf{r}) d\mathbf{r}

The usage pattern is as follows:

.. code-block:: python

    # Construct some potential, e.g. a hyperbolic well
    rsq = grid.points[:,0]**2 + grid.points[:,1]**2 + grid.points[:,2]**2
    pot = np.sqrt(1 + rsq)

    # Allocate an output array for the operator
    fock = lf.create_two_index()

    # Actual computation
    obasis.compute_grid_density_fock(grid.points, grid.weights, pot, fock)


Other chain rules are also implemented:

* :py:meth:`~horton.gbasis.cext.GOBasis.compute_grid_gradient_fock`

  .. math::

    \frac{\partial E[\nabla\rho]}{\partial D_{\nu\mu}} =
        \int \frac{\delta E[\nabla\rho]}{\delta \nabla\rho(\mathbf{r})}
        \cdot \left(
            \nabla b_\nu^*(\mathbf{r}) b_\mu(\mathbf{r}) +
            b_\nu^*(\mathbf{r}) \nabla b_\mu(\mathbf{r})
        \right) d\mathbf{r}

* :py:meth:`~horton.gbasis.cext.GOBasis.compute_grid_gradient_fock` just
  combines the density and gradient chain rules. This is more efficient than
  computing the separately.

* :py:meth:`~horton.gbasis.cext.GOBasis.compute_grid_kinetic_fock`

  .. math::

    \frac{\partial E[\tau]}{\partial D_{\nu\mu}} =
        \frac{1}{2}\int \frac{\delta E[\tau]}{\delta \tau(\mathbf{r})}
        \nabla b_\nu^*(\mathbf{r}) \nabla b_\mu(\mathbf{r}) d\mathbf{r}

  where :math:`\tau(\mathbf{r})` is the positive kinetic energy density:

  ..
    Mind the adjective "positive" in the following sentence. There are many
    choices for the kinetic energy density (which is essentially arbitrarily
    defined). Technically one should say "the nonnegative kinetic energy
    density" but I think most people call this the "positive" choice.
        ~ Paul W. Ayers

  .. math::

    \tau(\mathbf{r}) = \frac{1}{2} \sum_{\mu\nu} D_{\mu\nu} \nabla b_\nu^*(\mathbf{r}) \nabla b_\mu(\mathbf{r})

* :py:meth:`~horton.gbasis.cext.GOBasis.compute_grid_hessian_fock`

  .. math::

    \frac{\partial E[\nabla\nabla\rho]}{\partial D_{\nu\mu}} =
        \int \frac{\delta E[\nabla\nabla\rho]}{\delta \nabla\nabla\rho(\mathbf{r})}
        \colon \left(
            \nabla \nabla b_\nu^*(\mathbf{r}) b_\mu(\mathbf{r}) +
            2 \nabla b_\nu^*(\mathbf{r}) \nabla b_\mu(\mathbf{r}) +
            b_\nu^*(\mathbf{r}) \nabla \nabla b_\mu(\mathbf{r})
        \right) d\mathbf{r}

* :py:meth:`~horton.gbasis.cext.GOBasis.compute_grid_mgga_fock` just combines
  several of the previous chain rules: density, gradient, laplacian (trace of
  the Hessian) and kinetic energy density.
