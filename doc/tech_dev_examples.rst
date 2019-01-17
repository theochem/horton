..
    : HORTON: Helpful Open-source Research TOol for N-fermion systems.
    : Copyright (C) 2011-2019 The HORTON Development Team
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

Writing examples scripts
########################

The examples are located in the ``data/examples`` directory. The first level of
subdirectories is grouping the examples by category. Additional (small) input
files may be provided next to the example scripts. Make sure that each script
properly executes in a few seconds and that it does not require input that can't
be included in the examples directory (because it is too big or binary). Keep
the examples short and put excessive comments before each line to explain how
the example works. For each example, you should also add a corresponding test in
``horton/test/test_examples.py`` that just executes the example.
