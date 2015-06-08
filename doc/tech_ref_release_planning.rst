..
    : HORTON: Helpful Open-source Research TOol for N-fermion systems.
    : Copyright (C) 2011-2015 The HORTON Development Team
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

Release planning
################

Stable API and beta code
========================

In order to ensure that users and developers spend as little time as possible
performing maintenance work unrelated to research, the developers of HORTON have
committed to providing a stable Application Programmer Interface (API). This
means that most new releases of HORTON will remain compatible with your code and
will not change the behaviour of your code.

.. warning::

    *The stable API policy only applies to public members of code.* If you
    access private members of classes directly (i.e. _array in Matrix), you can
    expect breakage even between bugfix releases!

HORTON releases are given version numbers according to the `Semantic Versioning
Specification <http://semver.org/>`_, which directly relates to API stability:

- Major releases (i.e. 1.x.x -> 2.x.x) will not guarantee compatiblity with past
  versions.

  * We will wait a minimum of 18 months between major releases.
  * These releases will contain overarching restructuring of the code and
    renaming of functions and classes, including removal of deprecated code.


- Minor releases (i.e. 1.1.x -> 1.2.x) will remain compatible with past
  versions.

  * These releases will be used to introduce new features or tweaks to existing
    code.
  * New algorithms for existing features will be disabled by default in order to
    preserve old behaviour.


- Bugfix releases (i.e. 1.1.1 -> 1.1.2) will not make any changes to the code
  except to fix bugs.


- Occasionally, we will release features marked as *beta* to speed up sharing
  ideas.

  * They do not adhere to this stable API policy. We reserve the right to break
    API in beta features at any time.


Unstable code in HORTON 2.0.0
=============================

Some parts of the 2.0.0 version have no stable API yet because of known problems
or because the code was not fully reviewed yet.

* :py:mod:`horton.correaltedwfn`
* :py:mod:`horton.perturbation_theory`
* :py:mod:`horton.orbital_entanglement`
* :py:mod:`horton.orbital_utils`
* :py:mod:`horton.orbital_localization`

In future release, we will report in this sectoin the API changes for the
unstable parts.


Features in development
=======================

The following features are planned or are in a certain state of development.
If you are interested in testing or cobtributing in any of these areas, send a
mail to `the HORTON mailing list
<https://groups.google.com/forum/#!forum/horton-discuss>`_.

======== =======================================================================
 Target   Feature
======== =======================================================================
 2.1.0    MGGA Functionals
 2.1.0    Range-separated exchange integrals
======== =======================================================================
