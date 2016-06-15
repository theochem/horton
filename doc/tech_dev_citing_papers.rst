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

Citing scientific work
======================

Whenever you add a feature based on a scientific publication, it should be cited
properly:

1. Add an item to the file ``data/references.bib``. Include the ``doi`` if
   possible. (The ``url`` field can be used as an alternative if the ``doi`` is
   not available.) Maintain the chronological order. Do not use reference
   manager software to edit this file. Use a plain text editor instead and don't
   touch existing references if not needed.

2. Add ``log.cite('someref', 'a reason')`` to the code based on the publication, e.g.\
   ``log.cite('marques2012', 'using LibXC, the library of exchange and correlation functionals')``.
   This guarantees that paper is properly cited at the end of the HORTON screen
   output.
