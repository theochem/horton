.. _modphysham:

Model Hamiltonians
##################
.. contents::

At the moment Horton supports only the 1- dimensional (1-D) Hubbard model Hamiltonian.
This particular Hubbard model is the simplest model of interacting particles in a lattice,
with two terms in the Hamiltonian:

.. math::

    \hat{H}_{\rm Hub} = -t\sum_{j,\sigma} \left( a_{(j+1)\sigma}^{\dagger}a_{j\sigma}
    + a_{j\sigma}^{\dagger}a_{(j+1)\sigma} \right )
    +U\sum_j n_{j\uparrow} n_{j\downarrow},
where the first term (kinetic) represents nearest-neighbor hopping and the second term (electron repulsion)
is the repulsive on-site interaction. `t` and `U` are user specified parameters.

Doing calculations with model Hamiltonians one does not need to specify the molecule and obasis instances,
and starts directly by defining a ``LinalgFactory`` instance,

.. code-block:: python

   lf = DenseLinalgFactory(n)

where

   :lf: Only the ``DenseLinalgFactory`` is supported.
   :`n`: (int) is a number of sites

The number of doubly occupied sites is assigned by the ``AufbauOccModel`` model:

.. code-block:: python

   occ_model = AufbauOccModel(m)

where

   :lf: Only the ``DenseLinalgFactory`` is supported.
   :`m`: (int) is a number of doubly occupied sites


Horton allows you to use periodic boundary conditions (PBC) in the 1-D Hubbard model

.. code-block:: python

   modelham = Hubbard(pbc=True)

where

    :`pbc`: (logical) ``True`` is the default value.

The expansion coefficients are generated for the `n` sites

.. code-block:: python

    exp_alpha = lf.create_expansion(n)

where

   :`n`: (int) is a number of sites

It is necessary to define the overlap matrix:

.. code-block:: python

    olp = modelham.compute_overlap(lf)

where

   :lf: is an instance of  ``LinalgFactory``

The `t` parameter is embedded in the kinetic energy term:

.. code-block:: python

   kin = modelham.compute_kinetic(lf, -t)

where

   :lf: is an instance of  ``LinalgFactory``
   :`t`: (real) is the value of `t` parameter

The `U` parameter is assign to the electron repulsion energy term:

.. code-block:: python

    er = modelham.compute_er(lf, U)

where

   :lf: is an instance of  ``LinalgFactory``
   :`U`: (real) is the value of the `U` parameter

Finally, all terms of the 1-D Hubbard Hamiltonian are combined together:

.. code-block:: python

    terms = [
            RTwoIndexTerm(kin, 'kin'),
            RDirectTerm(er, 'hartree'),
            RExchangeTerm(er, 'x_hf'),
        ]

Example input file
===================

Restricted Hartree-Fock calculations using 1-D Hubbard model Hamiltonian with PBC
----------------------------------------------------------------------------------------------

In this particular example, the number of doubly-occupied sites is ``3``, the total number of sites is
``6``. The ``t`` parameter is set to -1, the ``U`` parameter is set to 2, and periodic boundary
conditions are employed.

.. code-block:: python

    from horton import *

    ###############################################################################
    ## Define Occupation model, expansion coefficients and overlap ################
    ###############################################################################
    lf = DenseLinalgFactory(6)
    occ_model = AufbauOccModel(3)
    modelham = Hubbard(pbc=True)
    exp_alpha = lf.create_expansion(6)
    olp = modelham.compute_overlap(lf)
    ###############################################################################
    # t-param, t = -1
    ###############################################################################
    kin = modelham.compute_kinetic(lf, -1)
    ###############################################################################
    # U-param, U = 2
    ###############################################################################
    er = modelham.compute_er(lf, 2)
    ###############################################################################
    ## Perform initial guess ######################################################
    ###############################################################################
    guess_core_hamiltonian(olp, kin, exp_alpha)
    terms = [
        RTwoIndexTerm(kin, 'kin'),
        RDirectTerm(er, 'hartree'),
        RExchangeTerm(er, 'x_hf'),
    ]
    ham = REffHam(terms)
    ###########################################################################################
    ## Do a Hartree-Fockk calculation #########################################################
    ###########################################################################################
    scf_solver = PlainSCFSolver()
    scf_solver(ham, lf, olp, occ_model, exp_alpha)
    energy = ham.compute()
