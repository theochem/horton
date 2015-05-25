.. _exportintegrals:

Dumping a Hamiltonian to a file
###############################


Export a Hamiltonian in ASCII format
====================================

Horton can export a Hamiltonian in ASCII format that uses the Molpro FCIDUMP file
convention (see :ref:`molprohamformat` for more details). If the Hamiltonian of
the system contains several one-electron terms (kinetic energy of electrons,
electron-nuclear attraction, etc.), they have to be combined to a single
contribution. For example, the kinetic energy and electron-nuclear
attraction terms can be combined as follows

.. code-block:: python

    one = kin.copy()
    one.iadd(na)

where ``kin`` and ``na`` are the electronic kinetic energy and the electron-nuclear
attraction terms, respectively. To write a Hamiltonian with the one- and two-electron integrals ``one`` and
``two`` to file, run :py:func:`horton.io.integrals_io.dump_fcidump`.

Note that you can specify an active space (it isn't used by default).
In case of an active space, the one-electron integrals :math:`t_{pq}` are

.. math::

    t_{pq} = \textrm{one}_{pq} + \sum_{i \in \textrm{ncore}} ( 2 \langle pi \vert qi \rangle - \langle pi \vert iq \rangle),

where :math:`\textrm{one}_{pq}` is the element :math:`pq` of :math:`\mathbf{one}` and
:math:`\langle pi \vert qi \rangle` is the appropriate two-electron integrals in physicist's notation.
The core energy of the active space is calculated as

.. math::

    e_{\rm core} = \textrm{ecore} + 2\sum_{i \in \textrm{ncore}} \textrm{one}_{ii} + \sum_{i, j \in \textrm{ncore}} (2 \langle ij \vert ij \rangle - \langle ij \vert ji \rangle)

where the two-electron integrals :math:`\langle pq \vert rs \rangle` contain only the
elements with active orbital indices :math:`p,q,r,s`. Note that only the symmetry-unique
elements of the one- and two-electron integrals are exported.


Example input files
===================

Exporting a Hamiltonian in ASCII format
---------------------------------------

In this example, we export the molecular Hamiltonian for the dinotrogen molecule
in the cc-pVDZ basis in the Hartree-Fock orbital basis. The Hamiltonian with all
active orbitals is written to ``FCIDUMP``, while the Hamiltonian in an active
space of CAS(8,8) is written to ``FCIDUMP8-8``.

.. code-block:: python

    #!/usr/bin/env python

    from horton import *

    ###############################################################################
    ## Set up molecule, define basis set ##########################################
    ###############################################################################
    mol = Molecule.from_file(context.get_fn('test/water.xyz'))
    obasis = get_gobasis(mol.coordinates, mol.numbers, 'cc-pvdz')
    ###########################################################################################
    ## Define Occupation model, expansion coefficients and overlap ############################
    ###########################################################################################
    lf = DenseLinalgFactory(obasis.nbasis)
    occ_model = AufbauOccModel(7)
    orb = lf.create_expansion(obasis.nbasis)
    olp = obasis.compute_overlap(lf)
    ###########################################################################################
    ## Construct Hamiltonian ##################################################################
    ###########################################################################################
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
    ###########################################################################################
    ## Perform initial guess ##################################################################
    ###########################################################################################
    guess_core_hamiltonian(olp, kin, na, orb)
    ###########################################################################################
    ## Do a Hartree-Fock calculation #########################################################
    ###########################################################################################
    scf_solver = PlainSCFSolver(1e-6)
    scf_solver(ham, lf, olp, occ_model, orb)
    ###########################################################################################
    ## Combine to single one-electron Hamiltonian #############################################
    ###########################################################################################
    one = kin.copy()
    one.iadd(na)

    ###########################################################################################
    ## Export Hamiltonian in Hartree-Fock molecular orbital basis (all orbitals active) #######
    ###########################################################################################
    dump_fcidump(lf, one, er, external['nn'], orb, 'FCIDUMP')

    ###########################################################################################
    ## Export Hamiltonian in Hartree-Fock molecular orbital basis for CAS(8,8) ################
    ###########################################################################################
    dump_fcidump(lf, one, er, external['nn'], orb, 'FCIDUMP8-8',
                      **{'nel': 8, 'ncore': 2, 'nactive': 8})
