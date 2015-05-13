.. _user_integration_grids:

How to work with integration grids in Horton
############################################

.. _user_integration_grids_specify:

Specifiying an integration grid
===============================

Horton primarily makes use of Becke-Lebedev grids for numerical integrals of
molecular volumes. If you are not familiar with this concept, then first read
Becke's paper: [becke1988_multicenter]_

In a nutshell, a Becke-Lebedev grid works as follows. Say, one is interested in
computing the integral of :math:`f(\mathbf{r})` over a molecular volume. In
practice, the integrand is often derived from the density and thus also contains
sharp spikes close to the atomic nuclei because. In order to integrate all these
unsmooth spikes properly, the molecular integral is first split into atomic
contributions:

.. math::
    \int f(\mathbf{r}) d\mathbf{r} = \sum_A \int w_A(\mathbf{r}) f(\mathbf{r}) d\mathbf{r}

where :math:`w_A(\mathbf{r})` is the atomic weight function for atom A. It is
1 close the nucleus of atom A and goes to zero inside the other atoms. Every
atomic integral is then computed on a grid in spherical coordinates. This is
typically a product grid, where different one-dimensional radial grids are
possible and the Lebedev-Laikov grids are always used for the angular part. In
the end, can always write the numerical integration as follows:

.. math::
    \int f(\mathbf{r}) d\mathbf{r} \approx \sum_{i=1}^{N_\text{grid}} w_i f(\mathbf{r}_i)

where :math:`N_\text{grid}` is the number of grid points, :math:`w_i` are the
integration grid weights and :math:`\mathbf{r}_i` are the integration grid
points.

Horton can automatically construct Becke-Lebedev integration grids for a given
molecular geometry. These are needed for a DFT computation or for an
atoms-in-molecules analysis. The default grid is constructed as follows:

.. code-block:: python

    grid = BeckeMolGrid(coordinates, numbers, pseudo_numbers)

where ``coordinates`` is an array with Cartesian coordinates of the atoms,
``numbers`` is an array with (integer) element numbers and ``pseudo_numbers``
is an array with floating point effective core charges. These arrays are always
available when one loads a molecule from a file. For example:

.. code-block:: python

    mol = Molecule.from_file('water.xyz')
    grid = BeckeMolGrid(mol.coordinates, mol.numbers, mol.pseudo_numbers)

A fourth argument can be provided to control the accuracy of the integration
grid, e.g.:

.. code-block:: python

    grid = BeckeMolGrid(coordinates, numbers, pseudo_numbers, 'fine')

The available levels of accuracy for the built-in integration grids are
documented here: :ref:`ref_grids`. One can also control in more detail the
radial and angular components of the integration grids. See the API
documentation of :py:class:`horton.grid.molgrid.BeckeMolGrid` for more details.


Computing a numerical integral involving the electron density
=============================================================

This section assumes that the following objects are already available:

* ``obasis``: an orbital basis set
* ``dm_full``: the spin-summed density matrix
* ``grid``: a Becke-Lebedev integration grid as introduce above.

If you are not familiar with the ``obasis`` and ``dm_full`` quantities, read
the sections :ref:`user_molecularham_basis` and :ref:`user_hf_dft`,
respectively. Note that the density matrix can also be loaded from a file
instead of computing it with Horton. This is demonstrated in the example at the
end of this section.

First one must evaluate the electron density on all grid points of the
integration grid. This can be done as follows:

.. code-block:: python

    rho = obasis.compute_grid_density_dm(dm_full, grid.points)

Several quantities can be evaluated on the grid, see
:py:meth:`~horton.gbasis.cext.GOBasis.compute_grid_density_dm`,
:py:meth:`~horton.gbasis.cext.GOBasis.compute_grid_gradient_dm`,
:py:meth:`~horton.gbasis.cext.GOBasis.compute_grid_kinetic_dm`,
:py:meth:`~horton.gbasis.cext.GOBasis.compute_grid_hartree_dm`,
:py:meth:`~horton.gbasis.cext.GOBasis.compute_grid_esp_dm`,
and
:py:meth:`~horton.gbasis.cext.GOBasis.compute_grid_orbitals_exp`
for more details. The total number of electrons can be easily computed as a
simple way to verify the accuracy of the grid:

.. code-block:: python

    print grid.integrate(rho)

Since ``rho`` is simply a Numpy array, it can be manipulated easily to compute
functions of the density, e.g.

.. code-block:: python

    print grid.integrate(rho**(4.0/3.0))

One can also use the array ``grid.points`` to evaluate expectation values
numerically, e.g. the following snippet evaluates the expectation value of
:math:`\vert\mathbf{r}\vert`:

.. code-block:: python

    r = (grid.points[:,0]**2 + grid.points[:,1]**2 + grid.points[:,2]**2)**0.5
    grid.integrate(rho, r)

As shown in the above snippet, the ``integrate`` method can take multiple
arguments that are all multiplied before integration.

The following script is a complete example that computes the expectation value
of :math:`\vert\mathbf{r}\vert` for an electron density loaded from a file.

.. literalinclude:: ../data/examples/grid/expectation_r.py
    :lines: 2-
    :caption: ../data/examples/grid/expectation_r.py


Constructing a one-body operator from a real-space potential
============================================================

This section assumes that the following objects are already available:

* ``obasis``: an orbital basis set
* ``lf``: an instance if ``DenseLinalgFactory`` or ``CholeskyLinalgFactory``
* ``dm_full``: the spin-summed density matrix
* ``grid``: a Becke-Lebedev integration grid as introduce above.

If you are not familiar with the ``obasis`` or ``lf``, go through the section
:ref:`user_molecularham_basis`. The density matrix can either be read from a
file or computed with Horton, see :ref:`user_hf_dft`.

Given a multplicative potential, one may write its expectation value as:

.. math::

    \langle V \rangle = \int \rho(\mathbf{r}) V(\mathbf{r}) d\mathbf{r}.

When the orbitals are expanded into a local basis set, this can be rewritten as:

.. math::

    \langle V \rangle = \sum_{\mu\nu} D_{\mu\nu} \mathcal{V}_{\nu\mu}

where :math:`D_{\mu\nu}` is the spin-summed density matrix. The matrix
:math:`\mathcal{V}_{\nu\mu}` is defined as

.. math::

    \mathcal{V}_{\nu\mu} = \int V(\mathbf{r}) b_\nu^*(\mathbf{r}) b_\mu(\mathbf{r}) d\mathbf{r}

where :math:`b_\mu(\mathbf{r})` are the orbital basis functions. Such matrices
can be constructed with the method
:py:meth:`horton.gbasis.cext.GOBasis.compute_grid_density_fock`. This is method
also useful when applying the chain rule to construct the contribution
to a Fock matrix from a density functional:

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


Similar methods,
:py:meth:`~horton.gbasis.cext.GOBasis.compute_grid_gradient_fock` and
:py:meth:`~horton.gbasis.cext.GOBasis.compute_grid_kinetic_fock`,
are available for the following two chain rules:

.. math::

    \frac{\partial E[\nabla\rho]}{\partial D_{\nu\mu}} & =
        \int \frac{\delta E[\nabla\rho]}{\delta \nabla\rho(\mathbf{r})}
        \cdot \left(
            \nabla b_\nu^*(\mathbf{r}) b_\mu(\mathbf{r}) +
            b_\nu^*(\mathbf{r}) \nabla b_\mu(\mathbf{r})
        \right) d\mathbf{r} \\
    \frac{\partial E[\tau]}{\partial D_{\nu\mu}} & =
        \frac{1}{2}\int \frac{\delta E[\tau]}{\delta \tau(\mathbf{r})}
        \nabla b_\nu^*(\mathbf{r}) b_\mu(\mathbf{r}) d\mathbf{r}

where :math:`\tau(\mathbf{r})` is the kinetic energy density:

.. math::

    \tau(\mathbf{r}) = \frac{1}{2} \sum_{\mu\nu} D_{\mu\nu} \nabla b_\nu^*(\mathbf{r}) \nabla b_\mu(\mathbf{r})
