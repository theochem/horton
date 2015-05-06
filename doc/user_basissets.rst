Basis sets
##########

.. contents::


Basis set format
==========================
Horton supports basis sets consisting of generally contracted Cartesian Gaussian functions up to a maximum angular momentum of six (h functions).
Horton is using the same basis set format as NWChem, and the basis sets can be downloaded from the EMSL webpage (https://bse.pnl.gov/bse/portal).


.. _horton-basis-library:

Basis set library
==========================
Most of the popular basis sets can be loaded automatically. A list of currently supported basis sets is tabulated below.

.. cssclass:: table-striped

======================== =========================== ===========================
Common basis name           Horton library name         Supported elements
======================== =========================== ===========================
STO-3G                     ``sto-3g``                         H--Kr
STO-6G                   ``sto-6g``                      H--Kr
3-21G                    ``3-21g``                       H--Kr
3-21G*                   ``3-21g*``                      Na--Ar
3-21+G*                  ``3-21+g*``                     Na--Ar
4-31G                    ``4-31g``                       H--Ne, P, S, Cl
6-31G                    ``6-31g``                       H--Zn
6-31G*                   ``6-31g*``                      H--Zn
6-31G**                  ``6-31g**``                     H--Zn
6-31+G                   ``6-31+g``                      H--Ca
6-31+G*                  ``6-31+g*``                     H--Ca
6-31++G*                 ``6-31++g*``                    H--Ca
6-31++G**                ``6-31++g**``                   H, Li--Ca
cc-pVDZ                  ``cc-pvdz``                     H--Ar, Ca-Kr
cc-pVTZ                  ``cc-pvtz``                     H--Ar, Ca-Kr
cc-pVQZ                  ``cc-pvqz``                     H--Ar, Ca-Kr
cc-pCVDZ                 ``cc-cpvdz``                    Li--Ar
cc-pCVTZ                 ``cc-cpvtz``                    Li--Ar
cc-pCVQZ                 ``cc-cpvqz``                    Li--Ar
aug-cc-pVDZ              ``aug-cc-pvdz``                 H--Ar, Sc--Kr
aug-cc-pVTZ              ``aug-cc-pvtz``                 H--Ar, Sc--Kr
aug-cc-pVQZ              ``aug-cc-pvqz``                 H--Ar, Sc--Kr
aug-cc-pCVDZ             ``aug-cc-cpvdz``                Li--Ar
aug-cc-pCVTZ             ``aug-cc-cpvtz``                Li--Ar
aug-cc-pCVQZ             ``aug-cc-cpvqz``                Li--Ar
def2-svpd                ``def2-svpd``                   H--Kr
def2-tzvp                ``def2-tzvp``                   H--Kr
def2-tzvpd               ``def2-tzvpd``                  H--Kr
def2-qzvp                ``def2-qzvp``                   H--Kr
def2-qzvpd               ``def2-qzvpd``                  H--Kr
======================== =========================== ===========================

These basis sets are taken from the EMSL library (https://bse.pnl.gov/bse/portal). When publishing results obtained with these basis sets,
please cite the following references [feller1996]_ and [Didier2007]_.

The basis sets are loaded using the following function call.

.. code-block:: python

   get_gobasis(coordinates, numbers, default, element_map, [index_map, pure])

with arguments

    :coordinates: (float)  A (N, 3) numpy array with Cartesian coordinates of the atoms (see :ref:`read-molgeometry`)
    :numbers: (int) A (N,) numpy vector with the atomic numbers (see :ref:`read-molgeometry`)
    :default: (str) The basis set name applied to each atom.

with optional arguments

    :element_map: (str) A dictionary with element names or numbers as keys, and basis sets as values. These specs override the default basis (default ``None``)
    :index_map: (str) A dictionary with atomic indexes (based on the order of the atoms) as keys and basis sets as values (default ``None``)
    :pure: (bol) By default pure (spherical) basis functions are used. Set this to false to switch to Cartesian basis functions (default ``True``)


Unique basis set for all atoms
==============================

Usually one wants to define an unique (the same) basis set for the whole system, this can be done by a function call.

.. code-block:: python

    obasis = get_gobasis(mol.coordinates, mol.numbers, 'cc-pvdz')

where `mol.coordinates` and `mol.numbers` are read from file (see :ref:`read-molgeometry`), and ``cc-pvdz`` is the cc-pVDZ basis set.


Specifying different basis sets for different atoms
===================================================

In some cases one wants to specify different basis sets for different atoms. For example, setting the 3-21G basis set for the hydrogen atom and the 6-31G basis set for the carbon atom, and STO-3G for all remainig atoms can be done as follows.

.. code-block:: python

    obasis = get_gobasis(mol.coordinates, mol.numbers, 'sto-3g',{'H':'3-21g','C':'6-31g'})

where `mol.coordinates` and `mol.numbers` are read from file (see :ref:`read-molgeometry`),  and ``sto-3g``, ``3-21g`` and ``6-31g`` are the basis set names (see :ref:`horton-basis-library`)

Alternatively, the same result can be obtained by substituting the H and C symbols with their atomic numbers

.. code-block:: python

    obasis = get_gobasis(mol.coordinates, mol.numbers, 'sto-3g',{'1':'3-21g','6':'6-31g'})
