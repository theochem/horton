Orbital entanglement analysis
#############################

.. contents::

.. _orbitalentanglement:

Orbital entanglement and orbital correlation
============================================

In quantum chemistry, the interaction of orbitals is commonly used to understand chemical processes. Examples are molecular orbital diagrams, frontier orbital theory, and ligand field theory. The concept of entanglement can be used to quantify the interaction of orbitals.

Specifically, the interaction of one orbital and all other orbitals can be measured by the single-orbital entropy :math:`s(1)_i` (or one-orbital entropy), which is calculated from the eigenvalues :math:`\omega_{\alpha,i}` of the one-orbital RDM

.. math::

    s(1)_i = -\sum_{\alpha=1}^4 \omega_{\alpha,i}\ln \omega_{\alpha,i}.

The one-orbital RDM is determined from an N-particle RDM by tracing out all other orbital-degrees of freedom except those of orbital *i*. This leads to a RDM whose dimension is equal to the dimension of the one-orbital Fock space (4-dimensional for spatial orbitals).

Similarly, the two-orbital entropy :math:`s(2)_{i,j}` quantifies the interaction of an orbital pair :math:`i,j` and all other orbitals. It is calculated from the eigenvalues of the two-orbital RDM :math:`\omega_{\alpha, i, j}` (with 16 possible states for spatial orbitals)

.. math::

    s(2)_{i,j} =-\sum_{\alpha=1}^{16} \omega_{\alpha, i, j} \ln \omega_{\alpha, i, j}.

The total amount of correlation between any pair of orbitals :math:`(i,j)` can be evaluated from the orbital-pair mutual information. The orbital-pair mutual information is calculated using the single- and two-orbital entropy and thus requires the one- and two-orbital RDMs,

.. math::

    I_{i|j} = \frac{1}{2} \big(s(2)_{i,j} - s(1)_{i} - s(1)_{j} \big) \big(1 - \delta_{ij}\big),

where :math:`\delta_{ij}` is the Kronecker delta.

Note that a correlated wavefunction is required to have non-zero orbital entanglement and correlation. In the case of an uncorrelated wavefunction (for instance, a single Slater determinant) the (orbital) entanglement entropy is zero.

For more information on orbital entanglement and correlation, see ref. [boguslawski2015a]_.


Supported features
==================

If not mentioned otherwise, the orbital entanglement module supports restricted orbitals, ``DenseLinalgFactory`` and ``CholeskyLinalgFactory``. The current version of Horton offers

1. Calculation of :ref:`single-orbital entropy and orbital-pair mutual information for a seniority-zero wavefunction <orbitalentanglementseniorityzero>`

Supported wavefunction models are

1. :ref:`AP1roG <introap1rog>` (seniority-zero wavefunction)


.. _orbitalentanglementseniorityzero:

Orbital entanglement and correlation for a seniority zero wavefunction
======================================================================

If you use this module, please cite [boguslawski2015a]_

How to set-up a calculation
---------------------------

To evaluate the single-orbital entropy and orbital-pair mutual information for a given wavefunction model, the corresponding 1- and 2-RDMs need to be calculated first. The 1- and 2-RDMs are used to evaluate the one- and two-orbital density matrices which eigenvalues are needed to calculate the single-orbital and two-orbital entropy.

Given the 1- and 2-RDMs, an instance of the ``OrbitalEntanglement`` class can be created,

.. code-block:: python

    entanglement = OrbitalEntanglementAp1rog(lf, one_dm, two_dm)

with arguments

    :lf: A linear algebra factory. One of ``DenseLinalgFactory``, ``CholeskyLinalgFactory``
    :one_dm: (``OneIndex`` instance) the 1-RDM
    :two_dm: (list of ``TwoIndex`` instances) the 2-RDM. The first list element contains the :math:`\Gamma_{pp}^{qq}` block, while the second list element contains the :math:`\Gamma_{pq}^{pq}` block.

To calculate the entanglement and correlation measures, use a function call,

.. code-block:: python

    entanglement()

By default, the single-orbital entropy and orbital-pair mutual information will be written to disk. The single-orbital entropy is stored in ``s1.dat``, while the orbital-pair mutual information is written to ``i12.dat``. ``s1.dat`` contains two columns, where the first column contains the orbital index and the second column the corresponding single-orbital entropy. ``i12.dat`` has 3 columns. The first two columns encode the orbital indices, while the third column contains the corresponding mutual information.

How to generate correlation diagrams
------------------------------------

Horton provides ``gnuplot`` scripts to generate the entanglement and correlation diagrams. All scripts are tested for ``gnuplot4.7``.

To generate the single-orbital entropy diagram, run

.. code-block:: bash

    build_so_entropy [init_index final_index]

where **init_index** and **final_index** are optional arguments. If provided, the single-orbital entropy will be plotted for orbital indices in the interval [init_index, final_index].

The orbital-pair mutual information plot can be generated by running

.. code-block:: bash

    build_mi cutoff [init_index final_index]

**cutoff** determines the lower cutoff value of the mutual information and must be given in orders of magnitude (1, 0.1, 0.01, 0.001, etc.). Orbital correlations that are smaller than **cutoff** will not be displayed in the mutual information diagram. As above, **init_index** and **final_index** are optional arguments. If provided, the mutual information will be plotted for orbital indices in the interval [init_index, final_index].


Example input files
===================

Orbital entanglement analysis of an AP1roG wavefunction
-------------------------------------------------------

This is a basic example on how to perform an orbital entanglement analysis in Horton. This script performs an orbital-optimized AP1roG calculation, followed by an orbital entanglement analysis of the AP1roG wavefunction for the water molecule using the cc-pVDZ basis set.

.. code-block:: python

    from horton import *
    ###############################################################################
    ## Set up molecule, define basis set ##########################################
    ###############################################################################
    mol = Molecule.from_file('mol.xyz')
    obasis = get_gobasis(mol.coordinates, mol.numbers, 'cc-pvdz')
    ###############################################################################
    ## Define Occupation model, expansion coefficients and overlap ################
    ###############################################################################
    lf = DenseLinalgFactory(obasis.nbasis)
    occ_model = AufbauOccModel(5)
    moceoff = lf.create_expansion(obasis.nbasis)
    olp = obasis.compute_overlap(lf)
    ###############################################################################
    ## Construct Hamiltonian ######################################################
    ###############################################################################
    kin = obasis.compute_kinetic(lf)
    na = obasis.compute_nuclear_attraction(mol.coordinates, mol.pseudo_numbers, lf)
    er = obasis.compute_electron_repulsion(lf)
    external = {'nn': compute_nucnuc(mol.coordinates, mol.pseudo_numbers)}
    terms = [
        RTwoIndexTerm(kin, 'kin'),
        RDirectTerm(er, 'hartree'),
        RExchangeTerm(er, 'x_hf'),
        RTwoIndexTerm(na, 'ne'),
    ]
    ham = REffHam(terms, external)
    ###############################################################################
    ## Perform initial guess ######################################################
    ###############################################################################
    guess_core_hamiltonian(olp, kin, na, moceoff)
    ###############################################################################
    ## Do a Hartree-Fock calculation ##############################################
    ###############################################################################
    scf_solver = PlainSCFSolver(1e-6)
    scf_solver(ham, lf, olp, occ_model, moceoff)
    ###############################################################################
    ## Combine one-electron integrals to single Hamiltonian #######################
    ###############################################################################
    one = kin.copy()
    one.iadd(na)

    ###############################################################################
    ## Do OO-AP1roG optimization ##################################################
    ###############################################################################
    ap1rog = RAp1rog(lf, occ_model)
    energy, g, l = ap1rog(one, er, external['nn'], moceoff, olp, True)

    ###############################################################################
    ## Calculate response density matrices ########################################
    ###############################################################################
    one_dm = lf.create_one_index()
    one_dm.assign(exp_alpha.occupations)
    twoppqq = lf.create_two_index()
    twopqpq = lf.create_two_index()
    twoppqq.compute_2dm_ap1rog(one_dm, g, l, 'ppqq')
    twopqpq.compute_2dm_ap1rog(one_dm, g, l, 'pqpq')

    ###############################################################################
    ## Do orbital entanglement analysis ###########################################
    ###############################################################################
    entanglement = OrbitalEntanglementAp1rog(lf, one_dm, [twoppqq,twopqpq])
    entanglement()
