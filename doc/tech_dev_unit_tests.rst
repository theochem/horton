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

Writing unit tests
##################

HORTON uses the `Nosetests <https://nose.readthedocs.org/en/latest/>`_
program to run all the unit tests. The goal of a unit test is to check whether
as small piece of code works as expected.


Running the tests
-----------------

The tests are run as follows (including preparation steps)::

    toony@poony ~/.../horton:master> ./cleanfiles.sh
    toony@poony ~/.../horton:master> ./setup.py build_ext -i
    toony@poony ~/.../horton:master> nosetests -v

This will run the tests with the version of HORTON in the source tree, i.e. not
the one that is installed with ``python setup.py install``. If you have not
removed any ``.py`` files, the script ``./cleanfiles.sh`` can be skipped. If you
did not change any low-level code, ``./setup.py build_ext -i`` can be skipped.

When working on a specific part of the code, it is often convenient to limit the
number of tests that are checked. The following runs only the tests in ``horton/test/test_cell.py``::

    toony@poony ~/.../horton:master> nosetests -v horton/test/test_cell.py

Within one file, you can also select one test function::

    toony@poony ~/.../horton:master> nosetests -v horton/test/test_cell.py:test_from_parameters3


Writing new tests
-----------------

All tests in HORTON are located in the directories ``horton/test`` and
``horton/*/test``. All module files containing tests have a filename that starts
with ``test_``. In these modules, all functions with a name that starts with
``test_`` are picked up by Nosetests. Tests that do not follow this convention
are simply ignored.

The basic structure of a test is as follows:

.. code-block:: python

    def test_sum():
        a = 1
        b = 2
        assert a+b == 3

HORTON currently contains many examples that can be used as a starting point
for new tests. The easiest way to write new tests is to just copy an existing
test (similar to what you have in mind) and start modifying it.

Most test packages in ``horton`` contain a ``common`` module that contains useful
functions that facilitate the development of tests. An important example is the
``check_delta`` function to test if analytical derivatives are properly
implemented. This is a simple example:

.. code-block:: python

    import numpy as np
    from horton.common import check_delta

    def test_quadratic():
        # a vector function that computes the squared norm divided by two
        def fun(x):
            return np.dot(x, x)/2

        # the gradient of that vector function
        def deriv(x):
            return x

        # the dimension of the vector x
        ndim = 5
        # the number of small displacements used in the test
        ndisp = 100
        # a reference point used for the test of the analytical derivatives
        x0 = np.random.uniform(-1, 1, ndim)
        # the small displacements, i.e. each row is one (small) relative vector
        dxs = np.random.uniform(1e-5, 1e5, (ndisp, ndim))

        check_delta(fun, deriv, x0, dxs)


Writing tests that need a temporary directory
---------------------------------------------

A context manager is implemented in ``horton.test.common`` to simplify tests
that need a temporary working directory. It can be used as follows:

.. code-block:: python

    from horton.test.common import tmpdir

    def test_something():
        with tmpdir('horton.somemodule.test.test_something') as dn:
            # dn is a string with the temporary directory name.
            # put here the part of the test that operates in the temporary directory.

On most systems, this temporary directory is a subdirectory of ``/tmp``. The
argument ``'horton.somemodule.test.test_something'`` will occur in the directory
name, such that it can be easily recognized if needed.


.. _tech_dev_unit_tests_random:

Writing tests that use random numbers
-------------------------------------

Tests that make use of random numbers can be problematic when they only fail sometimes for
very specific and rare values of the random numbers. To avoid issues with such corner
cases, one must fix the random seed in the tests as follows:

.. code-block:: python

    from horton.test.common import numpy_seed

    def test_foo():
        # Some deterministic test code here.
        # ...
        # The part of test that uses random numbers should be repeated several times
        # with different random numbers to make use of the randomness.
        for irep in xrange(100):
            # Fix the seed differently but determinstically at each repetition.
            with numpy_seed(irep):
                # Test code that uses numpy.random should go here.

The test ``test_tridiagsym_solve`` in ``horton/grid/test/test_cubic_spline.py`` is a
realistic example that properly uses random numbers.
