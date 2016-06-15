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

.. _user_ap1rog:

The AP1roG module
#################

Two-electron functions, called geminals, can be used to incorporate electron
correlation effects into the many-particle wavefunction. Horton supports a
special type of geminal-based wavefunction models, the antisymmetric product of
1-reference orbital geminals (AP1roG), which is equivalent to pair-coupled
cluster doubles. The AP1roG wavefunct ion effectively parameterizes the doubly
occupied configuration interaction wavefunction (DOCI), but requires only
mean-field computational cost in contrast to the factorial scaling of
traditional DOCI implementations [limacher2013]_. Currently the AP1roG module is
limited to closed-shell systems only.

.. _user_ap1rog_model:

The AP1roG model
================

The AP1roG wavefunction ansatz can be rewritten in terms of one-particle functions as a fully general pair-coupled-cluster wavefunction,

.. math::
    :label: ap1rog

    \vert \textrm{AP1roG}\rangle = \exp(\sum_{ia} c_i^a a^\dagger_a a^\dagger_{\bar{a}} a_{\bar{i}} a_i) \vert \Psi_0 \rangle

where :math:`a_p^{\dagger}`, :math:`a_{\bar{p}}^{\dagger}`, and :math:`a_p`, :math:`a_{\bar{p}}` are the electron creation and annihilation operators and :math:`p` and :math:`\bar{p}` denote :math:`\alpha` and :math:`\beta` spins, respectively. :math:`\vert \Psi_0 \rangle` is some independent-particle wave function (for instance, the Hartree−Fock determinant).
Indices :math:`i` and :math:`a` correspond to virtual and occupied orbitals with respect to :math:`\vert \Psi_0 \rangle`, :math:`P` and :math:`K` denote the number of electron pairs (:math:`P = N/2` with :math:`N` being the total number of electrons) and orbitals, respectively.
The geminal coefficient matrix (:math:`\mathbf{C}`) of AP1roG links the geminals with the underlying one-particle basis functions and has the following form,

.. math::
    :label: cia

    \mathbf{C}  =
    \begin{pmatrix}
      1      & 0       & \cdots & 0       & c_{1;P+1} & c_{1;P+2}&\cdots &c_{1;K}\\
      0      & 1       & \cdots & 0       & c_{2;P+1} & c_{2;P+2}&\cdots &c_{2;K}\\
      \vdots & \vdots  & \ddots & \vdots  & \vdots    & \vdots   &\ddots &\vdots\\
      0      & 0       & \cdots & 1       & c_{P;P+1} & c_{P;P+2}&\cdots & c_{P;K}
    \end{pmatrix}


The exponential form of eq. :eq:`ap1rog` assures the size extensivity of the geminal wavefunction, however, in order to ensure the size consistency, one has to optimize the orbitals (see [boguslawski2014a]_ and [boguslawski2014b]_). The simplest and most robust way is to use the variational orbital optimization (vOO-AP1roG) method implemented in HORTON (see :ref:`ooap1rog`).

Currently supported features
============================


If not mentioned otherwise, the AP1roG module supports spin-restricted orbitals and the ``DenseLinalgFactory`` and ``CholeskyLinalgFactory``. Specifically, the following features are provided:

    1. Optimization of AP1roG (eq. :eq:`ap1rog`) with (see :ref:`ooap1rog`) and without orbital optimization (see :ref:`ap1rog`) for a given Hamiltonian (see :ref:`preamble`).
    2. Variational orbital optimization and PS2c orbital optimization (see :ref:`keywords-oo-ap1rog` to choose the orbital optimizer)
    3. Calculation of response one- and two-particle reduced density matrices (see :ref:`responsedms`)
    4. Determination of AP1roG natural orbitals and occupation numbers (see :ref:`responsedms` and :ref:`natorb`)
    5. Calculation of the exact orbital Hessian (see :ref:`exacthessian`). Note that the orbital optimizer uses only a diagonal Hessian. The exact orbital Hessian can only be evaluated in combination with the ``DenseLinalgFactory``.

The AP1roG wave function and its response density matrices can then be used for post-processing. This version of HORTON offers:

    1. A posteriori addition of dynamic electron correlation using the perturbation module (see :ref:`pta` and :ref:`ptb` for documentation)
    2. Analysis of orbital correlations in the AP1roG wave function using the orbital entanglement module (see :ref:`orbital_entanglementseniorityzero` for documentation)
    3. Dump the Hamiltonian (collection of one- and two-electron integrals and the energy term due to core electrons and external potentials) in the AP1roG MO basis. The one- and two-electron integrals can be calculated for any pre-defined active space, that is, a selected number of active electrons and orbitals. The Hamiltonian is stored in the Molpro file format (see :ref:`hamiltonian_io` for documentation)


Input structure
===============

.. _preamble:

Getting started
---------------

To optimize an AP1roG wavefunction, the module requires a Hamiltonian and an initial guess for the orbitals (either an AO/MO coefficient matrix or an MO/MO coefficient matrix) as input arguments. HORTON provides different options for specifying the Hamiltonian and an orbital guess.

- The Hamiltonian is divided into three contributions: the one- and two-electron integrals as well as an external term (also referred to as core energy). Possible choices are:

    1. In-house calculation of the quantum chemical Hamiltonian expressed in the AO basis (kinetic energy of the electrons, electron-nuclear attraction, electron-electron repulsion, and nuclear-nuclear repulsion). All terms are calculated separately in HORTON (see :ref:`user_molecularham_matrix_elements` for documentation). Note, however, that all one-electron terms have to be combined into one single operator term. This can be done in the following way

        .. code-block:: python

            ###############################################################################
            ## Calculate kinetic energy (kin) and nuclear attraction (na) term ############
            ###############################################################################
            kin = obasis.compute_kinetic(lf)
            na = obasis.compute_nuclear_attraction(mol.coordinates, mol.pseudo_numbers, lf)
            ###############################################################################
            ## Combine one-electron integrals to single Hamiltonian #######################
            ###############################################################################
            one = kin.copy()
            one.iadd(na)


    2. In-house calculation of model Hamiltonians. Supported model Hamiltonians are summarized in :ref:`modphysham`. If the model Hamiltonian contains separate one-electron contributions, they have to be combined to a single operator as shown under point 1.


    3. External (one- and two-electron) integrals (in an orthonormal basis) and core energy can be read from file. The integral file must use the Molpro file format (see :ref:`hamiltonian_io` for more details). To load a Hamiltonian from file, run

        .. code-block:: python

             one, two, coreenergy = load_fcidump(lf, filename='./FCIDUMP')

      with arguments

        :lf: A linear algebra factory. Must be of type ``DenseLinalgFactory``. Note that ``CholeskyLinalgFactory`` is not supported

      and optional arguments

        :filename: (str) the filename of the fcidump file (default ``FCIDUMP``)

      The function ``load_fcidump`` has three return values; the one-electron integrals (``one``) stored as a ``TwoIndex`` object, the two-electron integrals (``two``) stored as a ``FourIndex`` object, and the core energy (``coreenergy``, float).

- A set of initial guess orbitals can be either generated in HORTON (including the AO overlap matrix) or read from disk (see :ref:`restart-ap1rog` to use orbitals generated in HORTON as initial guess). Examples for initial guess orbitals are:

    1. Restricted canonical Hartree-Fock orbitals (see :ref:`user_hf_dft`)

    2. Localized orbitals. HORTON supports Pipek-Mezey localization of canonical Hartree-Fock orbitals. See :ref:`localization` for documentation.

    3. If external integrals (expressed in an orthonormal basis) are used to define the Hamiltonian, the initial orbitals and the overlap matrix are the identity matrix and can be set as follows:

        .. code-block:: python

            orb = lf.create_expansion(nbasis)
            olp = lf.create_two_index(nbasis)
            olp.assign_diagonal(1.0)
            orb.assign(olp)

      where ``nbasis`` is the number of basis function (total number of orbitals in the active space).


.. _ooap1rog:

AP1roG with orbital optimization
--------------------------------

If you use this part of the module, please cite [boguslawski2014a]_ and [boguslawski2014b]_

.. _setup-oo-ap1rog:

How to set-up a calculation
^^^^^^^^^^^^^^^^^^^^^^^^^^^

After specifying a Hamiltonian and initial guess orbitals, you can create an instance of the ``RAp1rog`` class,

.. code-block:: python

    apr1og =  RAp1rog(lf, occ_model, npairs=None, nvirt=None)

with arguments

    :lf: A ``LinalgFactory`` instance
    :occ_model: (``AufbauOccModel`` instance) an Aufbau occupation model

and optional arguments

    :npairs: (int) number of electron pairs. If not specified, the number of pairs equals the number of occupied orbitals in the ``AufbauOccModel``
    :nvirt: (int) number of virtual orbitals. If not specified, the number of virtual orbitals is calculated from the total number of basis functions minus the number of electron pairs.

Note that no optional arguments need to be specified for the AP1roG model and the number of electron pairs and virtual orbitals is automatically determined using the ``AufbauOccModel``. A restricted, orbital-optimized AP1roG calculation can be initiated with a function call,

.. code-block:: python

    energy, c, l = ap1rog(one, two, external, orb, olp, scf, **keywords)

with arguments

    :one: (``TwoIndex`` instance) the one-electron integrals
    :two: (``FourIndex`` or ``CholeskyLinalgFactory`` instance) the two-electron integrals
    :external: (float) energy contribution due to an external potential, e.g., nuclear-nuclear repulsion term, etc.
    :orb: (``Expansion`` instance) the AO/MO or MO/MO coefficient matrix. It also contains information about orbital energies (not defined in the AP1roG model) and occupation numbers
    :olp: (``TwoIndex`` instance) the AO overlap matrix or, in case of an orthonormal basis, the identity matrix
    :scf: (boolean) if ``True``, orbitals are optimized

The keyword arguments are optional and contain optimization-specific options (like the number of orbital optimization steps, etc.) as well as orbital manipulation schemes (Givens rotations, swapping orbitals in reference determinant, etc.). Their default values are chosen to give reasonable performance and can be adjusted if convergence difficulties are encountered. All keyword arguments are summarized in the following section (:ref:`keywords-oo-ap1rog`).

The function call gives 3 return values,

    :energy: (float) the total AP1roG electronic energy (the **external** term included)
    :c: (``TwoIndex`` instance) the geminal coefficient matrix (without the diagonal occupied sub-block, see :ref:`user_ap1rog_model`)
    :l: (``TwoIndex`` instance) the Lagrange multipliers (can be used to calculated the response 1-RDM)

After the AP1roG calculation is finished (because AP1roG converged or the maximum number of iterations was reached), the orbitals and the overlap matrix are, by default, stored in a checkpoint file ``checkpoint.h5`` and can be used for a subsequent restart. Note that the geminal coefficient matrix and Lagrange multipliers are not stored after the calculation is completed.

.. _keywords-oo-ap1rog:

Summary of keyword arguments
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    :indextrans: (str) 4-index Transformation. Choice between ``tensordot`` (default) and ``einsum``. ``tensordot`` is faster than ``einsum``, but requires more memory. If ``DenseLinalgFactory`` is used, the memory requirement scales as :math:`2N^4` for ``einsum`` and :math:`3N^4` for ``tensordot``, respectively. Due to the storage of the two-electron integrals, the total amount of memory increases to :math:`3N^4` for ``einsum`` and :math:`4N^4` for ``tensordot``, respectively.

    :warning: (boolean) if ``True``, (scipy) solver-specific warnings are printed (default ``False``)

    :guess: (dictionary) initial guess specifications:

             :type: (str) guess type. One of ``random`` (random numbers, default), ``const`` (``1.0`` scaled by **factor**)
             :factor: (float) a scaling factor for the initial guess of type ``type`` (default ``-0.1``)
             :geminal: (1-dim np.array) external guess for geminal coefficients (default ``None``). If provided, **type** and **factor** are ignored. The elements of the geminal matrix of eq. :eq:`cia` have to be indexed in C-like order. Note that the identity block is not required. The size of the 1-dim np.array is thus equal to the number of unknowns, that is, :math:`n_{\rm pairs}*n_{\rm virtuals}`.
             :lagrange: (1-dim np.array) external guess for Lagrange multipliers (default ``None``). If provided, **type** and **factor** are ignored. The elements have to be indexed in C-like order. The size of the 1-dim np.array is equal to the number of unknowns, that is, :math:`n_{\rm pairs}*n_{\rm virtuals}`.

    :solver: (dictionary) scipy wavefunction/Lagrange solver:

             :wfn: (str) wavefunction solver (default ``krylov``)
             :lagrange: (str) Lagrange multiplier solver (default ``krylov``)

             Note that the exact Jacobian of **wfn** and **lagrange** is not supported. Thus, scipy solvers that need the exact Jacobian cannot be used. See `scipy root-solvers <http://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.optimize.root.html>`_ for more details.

    :maxiter: (dictionary) maximum number of iterations:

               :wfniter: (int) maximum number of iterations for the **wfn/lagrange** solver (default ``200``)
               :orbiter: (int) maximum number of orbital optimization steps (default ``100``)

    :thresh: (dictionary) optimization thresholds:

              :wfn: (float) optimization threshold for geminal coefficients and Lagrange multipliers (default ``1e-12``)
              :energy: (float) convergence threshold for energy (default ``1e-8``)
              :gradientnorm: (float) convergence threshold for norm of orbital gradient (default ``1e-4``)
              :gradientmax: (float) threshold for maximum absolute value of the orbital gradient (default ``5e-5``)

    :printoptions: (dictionary) print level:

              :geminal: (boolean) if True, geminal matrix is printed (default ``True``). Note that the identity block is omitted.
              :ci:  (float) threshold for CI coefficients (requires evaluation of a permanent). All coefficients (for a given excitation order) larger than **ci** are printed (default ``0.01``)
              :excitationlevel: (int) number of excited pairs w.r.t. the reference determinant for which the wavefunction amplitudes are reconstructed (default ``1``). At most, the coefficients corresponding to hextuply excited Slater determinants w.r.t the reference determinant can be calculated.

              Note that the reconstruction of the wavefunction amplitudes requires evaluating a permanent which is in general slow (in addition to the factorial number of determinants in the active space).

    :dumpci: (dictionary) dump Slater determinants and corresponding CI coefficients to file:

              :amplitudestofile: (boolean) write wavefunction amplitudes to file (default ``False``)
              :amplitudesfilename: (str) file name (default ``ap1rog_amplitudes.dat``)

    :stepsearch: (dictionary) optimizes an orbital rotation step:

              :method: (str) step search method used. One of ``trust-region`` (default), ``None``,  ``backtracking``
              :optimizer: (str) optimizes step to boundary of trust radius in ``trust-region``. One of ``pcg`` (preconditioned conjugate gradient), ``dogleg`` (Powell's single dogleg step), ``ddl`` (Powell's double-dogleg step) (default ``ddl``)
              :alpha: (float) scaling factor for Newton step. Used in ``backtracking`` and ``None`` method (default ``1.00``)
              :c1: (float) parameter used in the Armijo condition of ``backtracking`` (default ``1e-4``)
              :minalpha: (float) minimum step length used in ``backracking`` (default ``1e-6``). If step length falls below **minalpha**, the ``backtracking`` line search is terminated and the most recent step is accepted
              :maxiterouter: (int) maximum number of iterations to optimize orbital rotation step  (default ``10``)
              :maxiterinner: (int) maximum number of optimization steps in each step search (used only in ``pcg``, default ``500``)
              :maxeta: (float) upper bound for estimated vs. actual change in ``trust-region`` (default ``0.75``)
              :mineta: (float) lower bound for estimated vs. actual change in ``trust-region`` (default ``0.25``)
              :upscale: (float) scaling factor to increase trust radius in ``trust-region`` (default ``2.0``)
              :downscale: (float) scaling factor to decrease trust radius in ``trust-region`` (default ``0.25``)
              :trustradius: (float) initial trust radius (default ``0.75``)
              :maxtrustradius: (float) maximum trust radius (default ``0.75``)
              :threshold: (float) trust-region optimization threshold, only used in ``pcg`` (default ``1e-8``)

    :checkpoint: (int) frequency of checkpointing. If **checkpoint** > 0, orbitals (``orb.hdf5``) and overlap (``olp.hdf5``) are written to disk (default ``1``)

    :levelshift: (float) level shift of Hessian (default ``1e-8``). Absolute value of elements of the orbital Hessian smaller than **levelshift** are shifted by **levelshift**

    :absolute: (boolean), if ``True``, the absolute value of the orbital Hessian is taken (default ``False``)

    :sort: (boolean), if ``True``, orbitals are sorted according to their natural occupation numbers. This requires re-solving for the wavefunction after each orbital optimization step. Works only if **orbitaloptimizer** is set to ``variational`` (default ``True``)

    :swapa: (2-dim np.array) swap orbitals. Each row in **swapa** contains 2 orbital indices to be swapped (default ``np.array([[]])``)

    :givensrot: (2-dim np.array) rotate two orbitals using a Givens rotation. Each row in **givensrot** contains 2 orbital indices and the rotation angle in deg (default ``np.array([[]])``). Orbitals are rotated sequentially according to the rows in **givensrot**. If a sequence of Givens rotation is performed, note that the indices in **givensrot** refer to the already rotated orbital basis. If **givensrot** is combined with **swapa**, all orbitals are swapped prior to any Givens rotation

    :orbitaloptimizer: (str) switch between variational orbital optimization (``variational``) and PS2c orbital optimization (``ps2c``) (default ``variational``)


.. _restart-ap1rog:

How to restart
^^^^^^^^^^^^^^

To restart an AP1roG calculation (for instance, using the orbitals from a different molecular geometry as initial guess or from a previous calculation using the same molecular geometry), the molecular orbitals (of the previous wavefunction run) need to be read from disk,

.. code-block:: python

    old = IOData.from_file('checkpoint.h5')


If this checkpoint file is used as an initial guess for a new geometry, the
orbitals must be re-orthogonalized w.r.t. the new basis:

.. code-block:: python

    # It is assume that the following two variables available:
    #   olp: the overlap matrix of the new basis (geometry)
    #   orb: an instance of DenseExpansion to which the new orbitals will be
    #        written
    project_orbitals_ortho(old.olp, olp, old.exp_alpha, orb)

Then, generate an instance of the ``RAp1rog`` class and perform a function call:

.. code-block:: python

    apr1og =  RAp1rog(lf, occ_model)
    energy, c, l = ap1rog(one, two, external, orb, olp, True)

Note that all optional arguments have been omitted and ``orb`` was set to ``True``.


.. _responsedms:

Response density matrices
^^^^^^^^^^^^^^^^^^^^^^^^^

HORTON supports the calculation of the response 1- and 2-particle reduced density matrices (1-RDM and 2-RDM), :math:`\gamma_{pq}` and :math:`\Gamma_{pqrs}`, respectively. Since AP1roG is a product of natural geminals, the 1-RDM is diagonal and is calculated from

.. math::
    \gamma_p = \langle \Psi_0| (1+\hat{\Lambda}) a^\dagger_p a_p | \textrm{AP1roG} \rangle,

where :math:`\hat{\Lambda}` contains the deexcitation operator,

.. math::
    \hat{\Lambda} = \sum_{ia} \lambda_i^a (a^\dagger_i a^\dagger_{\bar{i}} a_{\bar{a}} a_a - c_i^a).

The response 1-RDM (a ``OneIndex`` instance) can be calculated in HORTON as follows

.. code-block:: python

    one_dm = lf.create_one_index()
    ap1rog.compute_1dm(one_dm, c, l, factor=2.0, response=True)

where ap1rog is an instance of the ``RAp1rog`` class and (see also :ref:`setup-oo-ap1rog` to get ``c`` and ``l``)
    :one_dm: (``OneIndex`` instance) output argument that contains the response 1-RDM in the _array attribute
    :c: (``TwoIndex`` instance) the geminal coefficient matrix (without the diagonal occupied sub-block, see :ref:`user_ap1rog_model`)
    :l: (``TwoIndex`` instance) the Lagrange multipliers

and optional arguments
    :factor: (float) a scaling factor for the 1-RDM. If ``factor=2.0``, the spin-summed 1-RDM is calculated, as :math:`\gamma_{p}=\gamma_{\bar{p}}` (default ``1.0``)
    :response: (boolean) if ``True``, the response 1-RDM is calculated (default ``True``)

The response 2-RDM is defined as

.. math::
    \Gamma_{pqrs} = \langle \Psi_0| (1+\hat{\Lambda})a^\dagger_p a^\dagger_{q}  a_{s} a_r| \textrm{AP1roG} \rangle.

In HORTON, only the non-zero elements of the response 2-RDM are calculated, which are :math:`\Gamma_{pqpq}=\Gamma_{p\bar{q}p\bar{q}}` and :math:`\Gamma_{p\bar{p}q\bar{q}}`. Specifically, the non-zero elements :math:`\Gamma_{pqpq}` and :math:`\Gamma_{ppqq}` (where we have omitted the information about electron spin) are calculated separately and stored as ``TwoIndex`` objects. Note that :math:`\gamma_p=\Gamma_{p\bar{p}p\bar{p}}`.

.. code-block:: python

    twoppqq = lf.create_two_index()
    twopqpq = lf.create_two_index()
    ###############################################################################
    ## Gamma_ppqq #################################################################
    ###############################################################################
    ap1rog.compute_2dm(twoppqq, one_dm, c, l, 'ppqq', response=True)
    ###############################################################################
    ## Gamma_pqpq #################################################################
    ###############################################################################
    ap1rog.compute_2dm(twopqpq, one_dm, c, l, 'pqpq', response=True)

with arguments (see again :ref:`setup-oo-ap1rog` to get ``c`` and ``l``)

    :twopqpq/twoppqq: (``TwoIndex`` instance) output argument that contains the response 2-RDM
    :one_dm: (``OneIndex`` instance) the response 1-RDM
    :c: (``TwoIndex`` instance) the geminal coefficient matrix (without the diagonal occupied sub-block, see :ref:`user_ap1rog_model`)
    :l: (``TwoIndex`` instance) the Lagrange multipliers

and optional arguments

    :response: (boolean) if ``True``, the response 1-RDM is calculated (default ``True``)

Note that, in HORTON, :math:`\Gamma_{p\bar{p}q\bar{q}} = 0 \, \forall \, p=q \in \textrm{occupied}` and :math:`\Gamma_{p\bar{q}p\bar{q}} =  0 \, \forall \, p=q \in \textrm{virtual}`.

.. _natorb:

Natural orbitals and occupation numbers
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If AP1roG converges, the final orbitals are the AP1roG natural orbitals and are stored in ``orb`` (see :ref:`ooap1rog` how to obtain ``orb``). The natural orbitals can be exported to the molden file format (see :ref:`ref_file_formats`) and visualized using, for instance, `Jmol <http://jmol.sourceforge.net>`_ or `VESTA <http://jp-minerals.org/vesta/en/>`_.

The natural occupation numbers, the eigenvalues of the response 1-RDM (see :ref:`responsedms` for how to calculate response RDMs) are stored in the ``occupations`` attribute (a 1-dim np.array) of ``orb`` and can be directly accessed after an AP1roG calculation using

.. code-block:: python

    orb.occupations

.. _exacthessian:

The exact orbital Hessian
^^^^^^^^^^^^^^^^^^^^^^^^^

Although the orbital optimizer uses a diagonal approximation to the exact orbital Hessian, the exact orbital Hessian can be evaluated after an AP1roG calculation. Note that this feature is only available for the ``DenseLinalgFactory``. The ``CholeskyLinalgFactory`` does not allow for the calculation of the exact orbital Hessian. Thus, this feature is limited by the memory bottleneck of the 4-index transformation of ``DenseLinalgFactory`` (see also ``indextrans`` in :ref:`keywords-oo-ap1rog`). To calculate the exact orbital Hessian, the one-electron (``one``) and two-electron integrals (``two``) need to be transformed into the AP1roG MO basis first,

.. code-block:: python

    onemo, twomo = transform_integrals(one, two, indextrans, orb)

with arguments

    :one: (``TwoIndex`` instance) the one-electron integrals
    :two: (``FourIndex`` instance) the two-electron integrals
    :indextrans: (str) the 4-index transformation, either ``tensordot`` (preferred) or ``einsum``
    :orb: (``Expansion`` instance) the AO/MO coefficient matrix

and return values

    :onemo: (list of ``TwoIndex`` instances, one element for each spin combination) the transformed one-electron integrals
    :twomo: (list of ``FourIndex`` instances, one element for each spin combination) the transformed two-electron integrals

This step can be skipped if the one- and two-electron integrals are already expressed in the (optimized) MO basis. The transformed one- and two-electron integrals (first element in each list) are passed as function arguments to the ``get_exact_hessian`` attribute function of ``RAp1rog`` which returns a 2-dim np.array with elements :math:`H_{pq,rs} = H_{p,q,r,s}`,

.. code-block:: python

    hessian = ap1rog.get_exact_hessian(onemo[0], twomo[0])

where ``ap1rog`` is an instance of ``RAp1rog`` (see :ref:`ooap1rog`). The exact orbital Hessian can be diagonalized using, for instance, the `np.linalg.eigvalsh routine <http://docs.scipy.org/doc/numpy/reference/generated/numpy.linalg.eigvalsh.html#numpy.linalg.eigvalsh>`_,

.. code-block:: python

    eigv = np.linalg.eigvalsh(hessian)

.. _ap1rog:

AP1roG without orbital optimization
-----------------------------------

If you use this part of the module, please cite [limacher2013]_

How to set-up a calculation
^^^^^^^^^^^^^^^^^^^^^^^^^^^

There are two different ways of performing a restricted AP1roG calculation without orbital optimization. Similar to the orbital-optimized counterpart, we have to create an instance of the ``RAp1rog`` class (using the same abbreviations for arguments as in :ref:`ooap1rog` and excluding all optional arguments):

.. code-block:: python

    apr1og =  RAp1rog(lf, occ_model)

The geminal coefficients can be optimized as follows

1. Use a function call and set ``orb=False``:

    .. code-block:: python

        energy, c = ap1rog(one, two, external, orb, olp, False, **keywords)

    where we have used the same abbreviations for function arguments as in :ref:`ooap1rog`. The keyword arguments contain optimisation-specific options and are summarized in :ref:`keywordsap1rog`.

    Note that the function call gives only 2 return values (the Lagrange multipliers are not calculated):

    :energy: (float) the total AP1roG electronic energy (the **external** term included)
    :c: (``TwoIndex`` instance) the geminal coefficient matrix (without the diagonal occupied sub-block, see :ref:`user_ap1rog_model`)

    In contrast to the orbital-optimized code, the orbitals and the overlap matrix are not stored to disk after the AP1roG calculation is finished.

2. Use the orbital-optimized version of AP1roG (see :ref:`ooap1rog`), but set ``orbiter`` to ``0``, which suppresses an orbital-rotation step:

    .. code-block:: python

        energy, c, l = ap1rog(one, two, external, orb, olp, True, **{
            'maxiter': {'orbiter': 0}
        })

    The function call gives 3 return values (total energy **energy**, geminal coefficients **c**, and Lagrange multipliers **l**). The Lagrange multipliers can be used for post-processing.

.. _keywordsap1rog:

Summary of keyword arguments
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    :indextrans: (str) 4-index Transformation. Choice between ``tensordot`` (default) and ``einsum``. ``tensordot`` is faster than ``einsum``, requires, however, more memory. If ``DenseLinalgFactory`` is used, the memory requirement scales as :math:`2N^4` for ``einsum`` and :math:`3N^4` for ``tensordot``, respectively. Due to the storage of the two-electron integrals, the total amount of memory increases to :math:`3N^4` for ``einsum`` and :math:`4N^4` for ``tensordot``, respectively.

    :warning: (boolean) if ``True``, (scipy) solver-specific warnings are printed (default ``False``)

    :guess: (dictionary) initial guess containing:

             :type: (str) guess type. One of ``random`` (random numbers, default), ``const`` (constant numbers)
             :factor: (float) a scaling factor for the initial guess of type ``type`` (default ``-0.1``)
             :geminal: (1-dim np.array) external guess for geminal coefficients (default ``None``). If provided, **type** and **factor** are ignored. The elements of the geminal matrix of eq. :eq:`cia` have to be indexed in C-like order. Note that the identity block is not required. The size of the 1-dim np.array is thus equal to the number of unknowns, that is, :math:`n_{\rm pairs}*n_{\rm virtuals}`.

    :solver: wfn/Lagrange solver (dictionary) containing:

             :wfn: (str) wavefunction solver (default ``krylov``)

             Note that the exact Jacobian of **wfn** is not supported. Thus, scipy solvers that need the exact Jacobian cannot be used. See `scipy root-solvers <http://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.optimize.root.html>`_ for more details.

    :maxiter: (dictionary) maximum number of iterations containing:

               :wfniter: (int) maximum number of iterations for the **wfn** solver (default ``200``)

    :thresh: (dictionary) optimization thresholds containing:

              :wfn: (float) optimization threshold for geminal coefficients (default ``1e-12``)

    :printoptions: print level; dictionary containing:

              :geminal: (boolean) if True, geminal matrix is printed (default ``True``). Note that the identity block is omitted.
              :ci:  (float) threshold for CI coefficients (requires evaluation of a permanent). All coefficients (for a given excitation order) larger than **ci** are printed (default ``0.01``)
              :excitationlevel: (int) number of excited pairs w.r.t. the reference determinant for which the wavefunction amplitudes are reconstructed (default ``1``). At most, the coefficients corresponding to hextuply excited Slater determinants w.r.t the reference determinant can be calculated.

              Note that the reconstruction of the wavefunction amplitudes requires evaluating a permanent which is in general slow (in addition to the factorial number of determinants in the active space).

    :dumpci: (dictionary) dump Slater determinants and corresponding CI coefficients to file:

              :amplitudestofile: (boolean) write wavefunction amplitudes to file (default ``False``)
              :amplitudesfilename: (str) file name (default ``ap1rog_amplitudes.dat``)

    :swapa: (2-dim np.array) swap orbitals. Each row in **swapa** contains 2 orbital indices to be swapped (default ``np.array([[]])``)

Troubleshooting in AP1roG-SCF calculations
==========================================

- **How to change the number of orbital optimization steps:**

  To increase the number of iterations in the orbital optimization, adjust the keyword ``maxiter`` (see :ref:`keywords-oo-ap1rog`):

  .. code-block:: python

      'maxiter': {'orbiter': int}

  where ``int`` is the desired number of iterations

- **The energy oscillates during orbital optimization:**

  The occupied-virtual separation breaks down and the reference determinant cannot be optimized. In some cases, fixing the reference determinant might accelerate convergence. However, the final solution might not be reasonable if the optimized geminal coefficient matrix contains elements that are significantly larger than 1.0 in absolute value. To fix the reference determinant, the ``sort`` keyword has to be set to ``False`` (see :ref:`keywords-oo-ap1rog`):

  .. code-block:: python

      'sort': False

- **The orbital optimization converges very, very slowly:**

  Usually, the orbital optimization converges fast around the equilibrium. For stretched distances (in the vicinity of dissociation, etc.) convergence can be very slow, especially if the final solution results in symmetry-broken orbitals. In such cases, the diagonal approximation to the Hessian is not optimal. However, the current version of HORTON does not support orbital optimization with the exact Hessian nor Hessian updates.

- **How to scan a potential energy surface**

  To accelerate convergence, restart from adjacent points on the potential energy surface (see :ref:`restart-ap1rog`). Using Hartree-Fock orbitals as initial guess might result in convergence difficulties and optimization problems of the reference determinant.

- **How to perturb the orbitals:**

  The initial guess orbitals can be perturbed using a sequence of Givens rotations (see also :ref:`keywords-oo-ap1rog`),

  .. code-block:: python

      'givensrot': np.array([[orbindex1a, orbindex1b, angle1],[orbindex2a, orbindex2b, angle2],...])

  where orbitals with indices **orbindex1a** and **orbindex1b** are rotated by angle **angle1**, etc. Givens rotations between orbital pairs can be used if, for instance, the orbital optimizer converges to a saddle point.


Example Python scripts
======================

Several complete examples can be found in the directory
``data/examples/ap1rog``. Three of these are discussed in the following
subsections.

The water molecule (a minimum input example)
--------------------------------------------

This is a basic example on how to perform an orbital-optimized AP1roG calculation in HORTON. This script performs an orbital-optimized AP1roG calculation on the water molecule using the cc-pVDZ basis set and RHF orbitals as initial orbitals.

.. literalinclude:: ../data/examples/ap1rog/water_minimal.py
    :caption: data/examples/ap1rog/water_minimal.py
    :lines: 2-


The water molecule (with all default keyword arguments)
-------------------------------------------------------

This is the same example as above, but all keyword arguments are mentioned explicitly using their default values.

.. literalinclude:: ../data/examples/ap1rog/water_default.py
    :caption: data/examples/ap1rog/water_default.py
    :lines: 2-


AP1roG with external integrals
------------------------------

This is a basic example on how to perform an orbital-optimized AP1roG calculation using one- and two-electron integrals from an external file. The number of doubly-occupied orbitals is ``5``, while the total number of basis functions is ``28``. See :ref:`modphysham`.


.. literalinclude:: ../data/examples/ap1rog/external_hamiltonian_n2_dense.py
    :caption: data/examples/ap1rog/external_hamiltonian_n2_dense.py
    :lines: 2-


AP1roG using model Hamiltonians
-------------------------------

This is a basic example on how to perform an orbital-optimized AP1roG calculation using 1-D Hubbard model
Hamiltonian. The number of doubly-occupied sites is ``3``, the total number of sites is ``6``. The ``t``
parameter is set to -1, the ``U`` parameter is set to 2, and periodic boundary conditions are employed.

.. literalinclude:: ../data/examples/ap1rog/hubbard.py
    :caption: data/examples/ap1rog/hubbard.py
    :lines: 2-

Note that for the Hubbard model, the external potential has to be set to ``0``,

.. code-block:: python

    energy, c, l = ap1rog(kin, two, 0, orb, olp, True)
