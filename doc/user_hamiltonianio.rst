How to export and import integrals
##################################

.. contents::

.. _readintegrals:

Importing integrals
===================


.. _exportintegrals:

Exporting integrals
===================


Example input files
===================

Importing a Hamiltonian
-----------------------

This is a basic example on how to read one- and two-electron integrals from an external file. The number of doubly-occupied orbitals is ``5``, while the total number of basis functions is ``28``.

.. code-block:: python

    from horton import *
    ###############################################################################
    ## Define number of occupied orbitals and total number of basis functions #####
    ###############################################################################
    nocc = 5
    nbasis = 28
    ###############################################################################
    ## Define Occupation model, expansion coefficients and overlap ################
    ###############################################################################
    lf = DenseLinalgFactory(nbasis)
    occ_model = AufbauOccModel(nocc)
    moceoff = lf.create_expansion(nbasis)
    olp = lf.create_two_index(nbasis)
    olp.assign_diagonal(1.0)
    moceoff.assign_diagonal(1.0)
    ###############################################################################
    ## Read Hamiltonian from file 'FCIDUMP' #######################################
    ###############################################################################
    one, two, core = read_integrals(lf, './FCIDUMP')


Exporting a Hamiltonian
-----------------------
