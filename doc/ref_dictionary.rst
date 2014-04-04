Horton dictionary of variables
##############################

For a full list of class names and variables, consult the :ref:`genindex`.

Rules
=====

The following naming rules are used in the Horton source code:

* Lower and upper case is used according to PEP8. In share: ``CamelCase`` for
  class names, ``lower_case`` for anything else. In lower-case names, underscore
  may be used to improve readability. For details, see:
  http://legacy.python.org/dev/peps/pep-0008/#prescriptive-naming-conventions.

* Except for local variables inside one function or method, variable names
  should at least have two characters.

* Lists or arrays are often written in the plural form, e.g. ``coordinates``.

* When a variable represents a number of elements of some kind, the prefix ``n``
  is used, e.g. ``natom``.

* For integer loop variables, use the prefix ``i``, e.g. ``iatom``. (Often a
  single ``i`` and ``j`` are still used for this purpose.)


Variable names and prefixes
===========================

:index:`at`
    Common prefix for a variable related to an atom. (Used when similar
    variables exist that are not related to atoms, e.g. ``atgrid``.)

:index:`coordintes`
    3D Cartesian coordinates of the nuclei. (Singular form is also used for the
    coordinate of a single atom.)

:index:`dm`
    A density matrix. (``dm_full``, ``dm_spin``, ``dm_alpha``, ``dm_beta`` for
    spin-summed, spin difference, alpha and beta density matrices.)

:index:`er` or `electron_repulsion`
    The electron repulsion two-body operator.

:index:`fock`
    A Fock matrix.

:index:`exp`
    An expansion of orbitals in a basis set, with orbital energies and
    occupation numbers. (``exp_alpha`` and ``exp_beta`` are typically used for
    alpha and beta orbitals, respectively.)

:index:`grid`
    The specification of an integration grid: the Cartesian coordinates of all
    the grid points and the corresponding integration weights.

:index:`ham`
    An (effective) hamiltonian.

:index:`kin` or :index:`kinetic`
    The kinetic energy operator.

:index:`lf`
    A ``LinalgFactory`` instance.

:index:`mol`
    A ``Molecule`` instance.

:index:`moldens`
    The spin-summed electron density (typically as an array of electron density
    values evaluated on a grid.)

:index:`na` or :index:`nuclear_attraction`
    The nuclear attraction operator.

:index:`numbers`
    An array with atomic numbers. (Singular form is also used for the
    atomic number of a single atom.)

:index:`log`
    The screen logger of Horton (See horton.log.)

:index:`obasis`
    An orbital basis.

:index:`occ_model`
    A model to assign occupation numbers to orbitals.

:index:`olp` or :index:`overlap`
    The overlap operator.

:index:`pseudo_numbers`
    Effective core charges. (Singular form is also used for the
    effective core charge of a single atom.)

:index:`scf_solver`
    An algorithm to optimize the orbitals as to minimize the energy of an
    effective Hamiltonian.

:index:`spindens`
    The alpha - beta electron density (typically as an array of electron density
    values evaluated on a grid.)

:index:`wpart` or :index:`cpart`
    A partitioning scheme.
