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


.. _tech_dev_checklist:

Branch review checklist
#######################

Before you do anything else, watch this: https://www.youtube.com/watch?v=wf-BqAjZb8M&index=2&list=PLdBBfnzuDrjFsiqkVw82l0VAOrUwGAsMn

This section is meant as a guideline for reviewing a pull request. However,
it is recommended that you also try to follow these guidelines while working on a new
branch, before you make a pull request.


Reviewing protocol
==================

Step 1. Automatic QA testing
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When a pull request is made, a bunch of QA tests are executed to detect common problems.
This reduces the burden of the reviewers but it is not a replacement for the remainder of
the reviewing protocol. We strive to avoid false-positives in the QA scripts such that any
error encountered must be fixed in the pull request. This implies that some tests (in
Pylint, Pycodestyle, ...) are disabled and must be inspected manually.

If the automatic QA tests do not pass, fix the issues with additional commits to the
branch of the pull request. These commits should appear automatically in the PR. When you
made whitespace errors or when you have to fix author names, rewrite the commits and push
them with the ``--force`` option.

If you want to run the QA tests locally, this can be done as follows:

1. Make sure all changes are committed in your feature branch and that the feature branch
   is checked out.

2. Optional: if you like, you can build several dependencies of HORTON, using exactly the
   same version as in the continuous integration (buildbot or Travis-CI):

   .. code-block:: bash

     ./tools/qa/deps/install_alldeps_twobranches.sh

2. Run the QA tests:

   .. code-block:: bash

     ./tools/qa/test_all_twobranches.sh

   This is the same script that is used in the continuous integration.

After the QA scripts pass on the build server, your code is ready to be reviewed by a human
being.


Step 2. Human review
~~~~~~~~~~~~~~~~~~~~

The sections below specify a minimum of mindful and mindless criteria that new code must
satisfy before it can be merged into the master branch. Some of these were checked by the
automatically but a large part still needs to be checked manually.

.. warning::

    "A foolish consistency is the hobgoblin of little minds." (PEP8)

    Humans tend to tackle the mindless criteria first, just because they are easier to
    check. However, when you overemphasize on the mindless ones, you'll be too busy to see
    the mindfull issues, even though these are still more important.

To the extent possible, all criteria that have to be checked have a label of the form
M?####. These make it easier for reviewers to comment on specific lines of code in PRs.

At least one person other than the author of the PR has to go through the change
set and check for conflicts with the criteria below. When issues are found, the author
should try to fix these with additional commits. In some cases, issues can only be
resolved by rewriting the history of the branch.


Step 3. Merging
~~~~~~~~~~~~~~~

Once the review is completed and all issues are properly solved, someone with write access
will check if the reviewing protocol was followed properly. If all is fine, the PR gets
merged into the master branch. PRs must be fast-forward merge-able. The author of the PR
is responsible for rebasing his/her branch if needed.


Mindful criteria
================

These are criteria that require the reviewer to actually understand the whole change-set.
They can't be formulated in terms of checking just a few lines at a time, which also means
such criteria are usually not testable with QA scripts.


Atomic commits
~~~~~~~~~~~~~~

* **MM101** Keep your Git commits as small as possible and avoid that unrelated changes
  are combined in one commit. This needs to be checked for all commits in a pull request.
  If needed, a history rewrite is in order.

* **MM102** When a commit fixes a bug, include a unit test in the same commit that would
  break without the fix.

* **MM103** When a commit adds a new feature, include a unit tests that validates the
  implementation.

More background can be found here: http://www.freshconsulting.com/atomic-commits/


Unit tests (other than coverage, see below)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* **MM201** The unit tests are divided into two categories: slow and fast ones, where fast
  means execution time below 1 second. The slow ones are marked with a descriptor in the
  unit tests source code, while all other ones are assumed to be fast. For example, a slow
  test is marked as follows:

  .. code-block:: python

      # In the beginning of the unit test module
      from nose.plugins.attrib import attr

      # The slow test
      @attr('slow')
      def test_something_slow():
          # slow test code

All new code must covered with fast unit tests, which is checked automatically. Also make
sure the fast unit tests follow these rules:

* **MM201** Code must be unit tested with all possible combinations of arguments and
  options or with a reasonable subset thereof.

* **MM202** Focus on testing the small components, rather than writing tests that depend
  on a lot of code.

* **MM203** Also write unit tests for simple functions.

* **MM204** Unit tests should not just run code but also check that the result is in line
  with expectations.

* **MM205** When tests use random numbers, use a fixed random seed, such that at every
  execution the same random numbers are used. This avoids that a test breaks every now and
  then. See :ref:`tech_dev_unit_tests_random`.

* **MM206** Avoid tests that just compare the output of a routine with the output from a
  previous initial version of that routine. Such tests need to be redone each time a bug
  is found in that routine.

* **MM207** Check implementations of analytic derivatives with the tools provided in
  HORTON. This is currently the ``check_delta`` function.

* **MM208** Use ``numpy.testing`` and ``nose.tools`` to make unit tests more readable.


Code structure
~~~~~~~~~~~~~~

* **MM301** Whenever using class inheritance, check if it would be more convenient to use
  composition instead.


Modularity
~~~~~~~~~~

Modular code is intuitively obtained through `separation of concerns
<https://en.wikipedia.org/wiki/Separation_of_concerns>`_. This means the following:

    **Different problems** are handled in different **maximally independent** pieces of
    code (modules).

There are two important parts in this rule:

* **Different problems.** You have to figure what this means in your case.
  Usually its straightforward, e.g. input/output code is different from actual
  computation code. You can (and should) obviously make the separation much more
  fine-grained. For example, different subproblems of a single computation can be
  separated.

* **Maximally independent.** Modules should depend on each other as little as possible,
  e.g.:

    * The public API of a module has to be as small as possible. If a module has a huge
      public API, it means that the user must have intimate knowledge of how the module
      works, which is non-modular.

    * The dependence graph of modules should not contain cycles. For example, when A
      depends on B and B depends on C, then C should not depend on A. It would be fine for
      A to depend on C.

    * Modules should have as little dependencies as possible.


Cython is evil
~~~~~~~~~~~~~~

Motivation:

* Cython is a constantly evolving language in which one can easily write very
  unreadable and dirty code. Good coding practices also evolve quickly as Cython is
  further developed.

* No QA tools (like Pycodestyle, Pydocstyle, Pylint, ...) exist that do some basic QA
  assurance. Everything has to be checked manually.

Work as follows:

* Either write Python or low-level code (C++, Fortran, ...)

* Only use Cython for wrapping the low-level code.


General mindless criteria
=========================

Miscellaneous (MG01##)
~~~~~~~~~~~~~~~~~~~~~~

* **MG0101** Never commit code that breaks unit tests This also means that all tests must
  pass on every commit. (TODO: this is currently not tested automatically. This would be
  doable with fast unit tests.)

* **MG0102** Never modify or tamper with QA scripts with no good reason.


Whitespace errors
~~~~~~~~~~~~~~~~~

The following are not allowed (and checked automatically):

* Lines ending with one or more whitespaces.
* Usage of tabs.
* Empty lines at the end of a file.


Comments (MG02##)
~~~~~~~~~~~~~~~~~

Python doc strings are great for documenting the API, but they are not sufficient to
document the internals of code. When someone wants to understand your code, additional
help is highly appreciated in the form of comments. At least include the following
information:

* **MG0201** A dictionary of local variables when variable names that do not speak
  for themselves. Such a dictionary can be written at the level of a single function or
  method, but it also makes sense at a class or even module level when certain variable
  names are consistently throughout. (Such consistency is preferable where possible.)

* **MG0202** When using variable names that correspond to symbols from equations in some
  paper, add a reference to that paper and the relevant equations. Also explicitly list
  the variables that correspond to symbols in the reference.

* **MG0203** Explain the role of groups of lines of source code in comments. Usually a
  single line does not need to be explained with a comment, unless it looks really
  cryptic. (Such cryptic statements should rather be avoided.) For example, it does not
  make sense to write the following comment because the Python code speaks for itself:

  .. code-block:: python

      # Select all strictly positive values from list l and assign the result to lpos.
      lpos = [value for value in l if value > 0]


Author name and e-mail (MG03##)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* **MG0301** Committer and Author e-mail addresses are checked automatically against the
  file AUTHORS in the source tree. This is just to make sure that everyone properly
  configures these settings in Git.

* **MG0302** Use names and e-mail addresses that you would use as corresponding author on a
  scientific paper.


Git commit message format (MG04##)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* **MG0401** Separate subject from body with a blank line
* **MG0402** Limit the subject line to 50 characters
* **MG0403** Capitalize the subject line
* **MG0404** Do not end the subject line with a period
* **MG0405** Use the imperative mood in the subject line
* **MG0406** Wrap the body at 72 characters
* **MG0407** Use the body to explain what and why vs. how

More background can be found here: http://chris.beams.io/posts/git-commit/

Consider setting the following environment variables such that ``vim`` is used as editor
for the commit messages. It offers syntax highlighting to facilitate writing good commit
messages:

.. code-block:: bash

    export VISUAL=vim
    export EDITOR="$VISUAL"


Units and unit conversion
~~~~~~~~~~~~~~~~~~~~~~~~~

TODO


Mindless Python criteria
========================

This section first goes over all the criteria that are checked automatically by different
programs, listing all enforced rules and exceptions. The last subsection discusses the
remaining criteria that have to be checked manually.

pycodestyle
~~~~~~~~~~~

See http://pycodestyle.readthedocs.io/

The complete list of error codes can be found here:
http://pycodestyle.readthedocs.io/en/latest/intro.html#error-codes

The following pycodestyle error codes are disabled.

* By default in pycodestyle (2.0.0):
    * **E121** (*) continuation line under-indented for hanging indent
    * **E123** (*) closing bracket does not match indentation of opening bracket’s line
    * **E126** (*) continuation line over-indented for hanging indent
    * **E226** (*) missing whitespace around arithmetic operator
    * **E241** (*) multiple spaces after ‘,’
    * **W503** (*) line break occurred before a binary operator

* Because they cause undesirable false positives:
    * **E127** continuation line over-indented for visual indent
    * **E128** continuation line under-indented for visual indent


Pydocstyle
~~~~~~~~~~

See http://pydocstyle.readthedocs.io/

All errors caught automatically by the Pydocstyle program must be fixed. Keep in mind that
this program does not cover all recommendations in PEP257.

A complete list of error messages can be found here:
http://pydocstyle.readthedocs.io/en/latest/error_codes.html

Pydocstyle is executed with the default settings, except that the following is disabled:

    * **D103**: Missing docstring in public function

This is already checked by Pylint and is not enforced for ``test_*`` functions.
Pycodestyle cannot yet be configured to ignore test functions. (It can only ignore test
files which is not fine-grained enough for our purposes.)


PyLint
~~~~~~

See https://www.pylint.org/.

The complete list of error messages can be found here:
https://docs.pylint.org/features.html

The following messages are excluded by default: I0020, I0021, W0704. (It is not clear what
these stand for. They are not documented in Pylint.)

The following messages are excluded by default but are activated in our case (related to Python
3): E1601, E1602, E1603, E1604, E1605, E1606, E1607, E1608, W1601, W1602, W1603, W1604,
W1605, W1606, W1607, W1608, W1609, W1610, W1611, W1612, W1613, W1614, W1615, W1616, W1617,
W1618, W1619, W1620, W1621, W1622, W1623, W1624, W1625, W1626, W1627, W1628, W1629, W1630,
W1632, W1633, W1634, W1635, W1636, W1637, W1638, W1639, W1640

The following are excluded because we don't consider them the be fatal:

* **C0103**: invalid-name. Invalid %s name “%s”%s Used when the name doesn't match
  the regular expression associated to its type (constant, variable, class...).
* **I0011**: locally-disabled. Used when an inline option disables a message or a messages
  category.
* **W0613**: unused-argument. Unused argument %r Used when a function or method argument
  is not used.

The following is disabled to allow access to the protected members of the
(not-so-well-designed) Matrix classes:

* **W0212**: protected-access. Access to a protected member %s of a client class Used when
  a protected member (i.e. class member with a name beginning with an underscore) is
  access outside the class or a descendant of the class where it’s defined.

The following are excluded due false positives:

* **E0611**: no-name-in-module. No name %r in module %r Used when a name cannot be found
  in a module.
* **E1136**: Value ‘%s’ is unsubscriptable emitted when a subscripted value
  doesn’t support subscription(i.e. doesn’t define __getitem__ method)
* **E1101**: no-member. %s %r has no %r member Used when a variable is accessed for an
  unexistent member.
* **R0201**: no-self-use. Method could be a function Used when a method doesn't use its
  bound instance, and so could be written as a function.
* **C0411**: wrong-import-order. %s comes before %s Used when PEP8 import order is not
  respected (standard imports first, then third-party libraries, then local imports)
* **W0621**: Redefining name %r from outer scope (line %s) Used when a variable’s
  name hide a name defined in the outer scope.

The PyLint settings used by the QA scripts can be found in ``tools/qa/pylintrc``. Some
of the non-default settings in that file include:

* No doc strings are required for unit tests, i.e. functions starting with ``test_``.

* Variables that are intentionally unused should get the prefix ``_``. These are dummy
  variables. These may be useful when receiving return values that are not used, e.g.

  .. code-block:: python

    a, _b, c = some_function()  # No intent to use _b

    for dirpath, _dirnames, filenames in os.walk(source_directory):
        # No intent to use _dirnames in the body of the for loop

* The parameters for the design checks are significantly relaxed.


Code coverage by (fast) unit tests
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

QA scripts will check if new Python code is touched by unit tests.


Mindless Python criteria to be checked manually
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The following sources were used to compile the list of criteria below, sometimes making
verbatim copies:

* PEP8: http://www.python.org/dev/peps/pep-0008/ (version 01-Aug-2013)
* PEP257: http://www.python.org/dev/peps/pep-0257/ (vesion 13-Jun-2001)
* GPSG: https://google.github.io/styleguide/pyguide.html (Revision 2.59)

For every criterion below, the source is mentioned. If the source is prefixed with a
twidle, it means that we intentionally deviate from recommendations given in the source.
When no source is mentioned, the criteria are specific to HORTON.

* **MP01##** Docstrings

    * **MP0101** (~PEP257) Usage information of a script does not have to be listed in its
      module docstring. We use ``argparse`` instead, which also produces nice usage
      documentation when the script is called with ``-h``.

    * **MP0102** (~PEP257) A list of classes, functions, etc in the module docstring is not
      required as such tables of content are generated automatically by Sphinx.

    * **MP0103** (PEP257) Module docstrings should start with a short title followed by an
      empty line. (This is also assumed by scripts that generate the API reference
      documentation.)

    * **MP0104** (PEP257) A Module docstring should give some basic background on the
      module and include some example usage.

    * **MP0105** (PEP257) Class doc strings should explain the purpose and behavior of a
      class.

    * **MP0106** (PEP257) A base class docstring must explain how to implement derived
      classes.

    * **MP0107** (PEP257) Function and method docstrings must use an imperative mood in
      their first line.

    * **MP0108** Docstrings must be written in `RestructedText
      <http://sphinx-doc.org/rest.html>`_.

    * **MP0109** `Numpy docstring conventions
      <https://github.com/numpy/numpy/blob/master/doc/HOWTO_DOCUMENT.rst.txt>`_ must be
      followed.


* **MP02##** Import conventions

    * **MP0201** NumPy, H5Py, PyPlot and SciPy packages must be imported in the following
      way:

      .. code-block:: python

            # Must be on separate lines
            import numpy as np
            import h5py as h5
            import matplotlib.pyplot as pt
            # scipy subpackages have to import separately
            from scipy import whatever

    * **MP0202** (PEP8) Wildcard imports are (only) allowed in two situations:

        1. In ``__init__.py`` files to republish the API of submodules and subpackages.

        2. In unit tests, one may write:

           .. code-block:: python

                from horton import *

           This tests if everything in HORTON can be properly imported.

    * **MP0203** (PEP8) Put any relevant ``__all__`` specification directly after the
      imports.

    * **MP0204** (beyond PEP8) Never use relative imports.

    * **MP0205** (PEP8) When importing a class from a class-containing module, it's usually okay
      to spell this:

      .. code-block:: python

        from myclass import MyClass
        from foo.bar.yourclass import YourClass

      If this spelling causes local name clashes, then spell them

      .. code-block:: python

        import myclass
        import foo.bar.yourclass

      and use ``myclass.MyClass`` and ``foo.bar.yourclass.YourClass``.


* **MP03##** Naming conventions

    * **MP0301** Use self-explaining variable, function, method, class, module names where
      reasonable.

    * **MP0302** For integer quantities of something: ``nsomething``, e.g. ``natom``.

    * **MP0303** Loop variable in loop over a number of things: ``isomething``, e.g. for
      number of atoms:

      .. code-block:: python

          for iatom in xrange(natom):
              # loop content

    * **MP0304** Use plural for arrays, e.g. ``coordinates``

    * **MP0305** When using ``*`` and ``**`` constructs in python to allow for an arbitrary
      number of (keyword) arguments to functions, then always use the names ``*args``
      and/or ``**kwargs``.

    * **MP0306** (PEP8) Use ``_single_leading_underscore`` as a weak "internal use"
      indicator. E.g. ``from M import *`` does not import objects whose name starts with
      an underscore.

    * **MP0307** (PEP8) Use ``single_trailing_underscore_`` as to avoid conflicts with
      Python keyword.

    * **MP0308** (PEP8) ``__double_leading_and_trailing_underscore__`` are "magic" objects
      or attributes that live in user-controlled namespaces. E.g. __init__ , __import__ or
      __file__ . Never invent such names; only use them as documented.

    * **MP0309** (PEP8) Modules should have short, all-lowercase names. Underscores can
      be used in the module name if it improves readability. Python packages should also
      have short, all-lowercase names, although the use of underscores is discouraged.

    * **MP0310** (PEP8) Class names should normally use the ``CapWords`` convention.
      When abbreviations and acronyms are used in a name, capitalize them, e.g.
      ``CPPCheck``.

    * **MP0311** (PEP8) Use the suffix "Error" on your exception names (if the exception
      actually is an error).

    * **MP0312** (PEP8) Function, method and variable names should be lowercase, with words
      separated by underscores as necessary to improve readability.

    * **MP0313** (PEP8) Always use self for the first argument to instance methods.

    * **MP0314** (PEP8) Always use cls for the first argument to class methods.

    * **MP0315** (PEP8) Constants are usually defined on a module level and written in all
      capital letters with underscores separating words. Examples include ``MAX_OVERFLOW``
      and ``TOTAL``.


* **MP04##** Cite papers where appropriate. Whenever you add a feature based on a
  scientific publication, it should be cited properly:

    * **MP0401** Add an item to the file data/references.bib. Use a lowercase bibtex key,
      following the lastnameyear convention. Include the doi if possible.
      (The url field can be used as an alternative if the doi is not available.)
      Maintain chronological order and alphabetical order within one year.

    * **MP0402** Cite the references in the HORTON documentation as follows:
      ``[lastnameyear]_``

    * **MP0403** Add ``log.cite('someref', 'a reason')`` to the code based on the
      publication, e.g. ``log.cite('marques2012', 'using LibXC, the library of exchange
      and correlation functionals').``

    * **MP0404** No begging for citations in the output or documentation. Go beyond
      self-citations. Try to be informative and neutral.


* **MP05##** Code formatting (wrapping, indentation, whitespace, ...)

    * **MP0501** (PEP8) Smart wrapping with parenthesis is preferred over wrapping with a
      backslash. However, in some cases, this may not be appropriate, like in ``if``,
      ``with``, ``assert`` statements etc.

    * **MP0502** (PEP8) Extra blank lines may be used (sparingly) to separate groups of
      related functions. Blank lines may be omitted between a bunch of related one-liners
      (e.g. a set of dummy implementations).

    * **MP0503** (PEP8) Use blank lines in functions, sparingly, to indicate logical
      sections.

    * **MP0504** (PEP8) The first or second line of each Python file must contain
      ``# -*- coding: UTF-8 -*-"``. (Only needed for Python 2.)

    * **MP0505** (PEP8) If operators with different priorities are used, consider adding
      whitespace around the operators with the lowest priority(ies). Use your own
      judgment; however, never use more than one space, and always have the same amount of
      whitespace on both sides of a binary operator.

      Yes:

      .. code-block:: python

        i = i + 1
        submitted += 1
        x = x*2 - 1
        hypot2 = x*x + y*y
        c = (a+b) * (a-b)

      No:

      .. code-block:: python

        i=i+1
        submitted +=1
        x = x * 2 - 1
        hypot2 = x * x + y * y
        c = (a + b) * (a - b)

    * **MP0506** (PEP8) Don't use spaces around the ``=`` sign when used to indicate a
      keyword argument or a default parameter value.

      Yes:

      .. code-block:: python

        def munge(input: AnyStr):
        def munge(sep: AnyStr = None):
        def munge() -> AnyStr:
        def munge(input: AnyStr, sep: AnyStr = None, limit=1000):

      No:

      .. code-block:: python

        def munge(input: AnyStr=None):
        def munge(input:AnyStr):
        def munge(input: AnyStr)->PosInt:

    * **MP0507** (PEP8) Do use spaces around the = sign of an annotated function definition.
      Additionally, use a single space after the : , as well as a single space on either
      side of the -> sign representing an annotated return value.

    * **MP0508** (GPSG) Every file should contain license boilerplate.

    * **MP0509** Scripts (and scripts only) should have ``#!/usr/bin/env python`` as their
      first line.


* **MP06##** Follow the PEP8 rules for comments, except for Strunk and White if you know
  better: https://www.python.org/dev/peps/pep-0008/#comments. (See also block comments and
  inline comments.) TODO: make list of bullet points instead of just linking to PEP8


* **MP07##** API

    * **MP0701** (PEP8, GPSG) Do not define public ``get_*`` or ``set_*`` methods in a
      class that involve litte computation. Make these methods non-public (prefix with
      underscore) and wrap them in a property instead. If these methods involve
      significant computation, they are fine as a method but try to find a better name

    * **MP0702** Follow the PEP8 rules given here:
      https://www.python.org/dev/peps/pep-0008/#public-and-internal-interfaces
      TODO: make list of bullet points instead of just linking to PEP8

    * **MP0703** (PEP8) All of the following applies, except that we completely dissalow
      name mangling: https://www.python.org/dev/peps/pep-0008/#designing-for-inheritance
      TODO: make list of bullet points instead of just linking to PEP8


* **MP08##** Exception handling

    * **MP0801** (PEP8) Derive exceptions from ``Exception`` rather than ``BaseException``.

    * **MP0802** (PEP8) When raising an exception in Python 2, use raise
      ``ValueError('message')`` instead of the older form raise ``ValueError, 'message'``.
      The latter form is not legal Python 3 syntax.

    * **MP0803** (PEP8) When catching exceptions, mention specific exceptions whenever
      possible instead of using a bare ``except:`` clause.

    * **MP0804** (PEP8) When binding caught exceptions to a name, prefer the explicit name
      binding syntax added in Python 2.6:

      .. code-block:: python

        try:
            process_data()
        except Exception as exc:
            raise DataProcessingFailedError(str(exc))

      This is the only syntax supported in Python 3, and avoids the ambiguity problems
      associated with the older comma-based syntax.

    * **MP0805** (PEP8) Additionally, for all try/except clauses, limit the try clause to
      the absolute minimum amount of code necessary. This avoids masking bugs.

      Yes:

      .. code-block:: python

        try:
            value = collection[key]
        except KeyError:
            return key_not_found(key)
        else:
            return handle_value(value)

      No:

      .. code-block:: python

        try:
            # Too broad!
            return handle_value(collection[key])
        except KeyError:
            # Will also catch KeyError raised by handle_value()
            return key_not_found(key)

    * **MP0806** (GPSG) When the built-in exceptions seem inappropriate or too vague (e.g.
      like ``RuntimeError``), modules should defined their own exceptions. These
      exceptions must be defined in the module where they are used.


* **MP09##** Boolean expressions, implicit True/False

    * **MP0905** (PEP8) Don't compare boolean values to True or False using == .

      .. code-block:: python

        # Yes:
        if greeting:
        # No:
        if greeting == True:
        # Worse:
        if greeting is True:

    * **MP0906** (~PEP8, ~GPSG) For sequences, (strings, lists, tuples), DO NOT use the
      fact that empty sequences are false.

        * This is recommended in PEP8 but it is not very readable as one has to know the
          type of ``seq`` to figure out what is going on:

          .. code-block:: python

            if not seq:
            if seq:

        * Not good either:

          .. code-block:: python

            if len(seq):
            if not len(seq):

        * Recommended explicit form:

          .. code-block:: python

            if len(seq) > 0:
            if len(seq) == 0:

    * **MP0906** (GPSG) When handling integers, implicit false may involve more risk than
      benefit (i.e., accidentally handling None as 0).


* **MP99##** Miscellaneous

    * **MP9901** (PEP8) When a resource is local to a particular section of code, use a
      ``with`` statement to ensure it is cleaned up promptly and reliably after use.

    * **MP9902** (PEP8) Be consistent in return statements. Either all return statements
      in a function should return an expression, or none of them should. If any return
      statement returns an expression, any return statements where no value is returned
      should explicitly state this as ``return None``, and an explicit return statement
      should be present at the end of the function (if reachable).

    * **MP9903** (PEP8) Use string methods instead of the ``string`` module.

    * **MP9904** (PEP8) Use ``.startswith()`` and ``.endswith()`` instead of string
      slicing to check for prefixes or suffixes.

    * **MP9905** (GPSG) Use default iterators and operators for types that support them, like
      lists, dictionaries, and files. The built-in types define iterator methods, too.
      Prefer these methods to methods that return lists, except that you should not mutate
      a container while iterating over it.

      Yes:

      .. code-block:: python

          for key in adict: ...
          if key not in adict: ...
          if obj in alist: ...
          for line in afile: ...
          for k, v in dict.iteritems(): ...

      No:

      .. code-block:: python

          for key in adict.keys(): ...
          if not adict.has_key(key): ...
          for line in afile.readlines(): ...

    * **MP9906** (GPSG) Avoid unreadable functional programming constructs, e.g. because
      they do not fit on one or two lines. This includes list comprehensions, lambda
      functions and inline conditionals.

    * **MP9907** (GPSG, PEP8) If a class inherits from no other base classes, explicitly
      inherit from object. This also applies to nested classes.

    * **MP9908** (GPSG) Your code should always check if ``__name__ == '__main__'`` before
      executing your main program so that the main program is not executed when the module
      is imported. Just call the ``main`` function in this if clause instead of adding a
      lot of code in the global scope.

      .. code-block:: python

        def main():
            # ...

        if __name__ == '__main__':
            main()


Mindless C++ criteria
=====================

All C++ code should make use of the C++11 standard. Automatic checks are only applied to
manually written C++ code. Autogenerated code is excluded from such tests.


CPPCheck
~~~~~~~~

See http://cppcheck.sourceforge.net/

CPPCheck is executed with all checks enabled and with the C++11 flag. The following
exceptions are added due to false positives:

- ``missingIncludeSystem``
- ``unusedFunction``

All other errors must be fixed.


CPPLint
~~~~~~~

See https://github.com/google/styleguide/tree/gh-pages/cpplint

CPPLint is executed with the default settings, except that the maximum line length is set
to 100. All errors must be fixed.


Manual checks
~~~~~~~~~~~~~

The following points should be checked manually. These are taken from the Google C++ Style
Guide (GCSG). See https://google.github.io/styleguide/cppguide.html

* **MC00##** Header files

    * **MC0001** (GCSG) `#define guard <https://google.github.io/styleguide/cppguide.html#The__define_Guard>`_
      All header files should have #define guards to prevent multiple inclusion. The
      format of the symbol name should be ``<PROJECT>_<PATH>_<FILE>_H_``.

    * **MC0004** (GCSG) `Names and Order of Includes <https://google.github.io/styleguide/cppguide.html#Names_and_Order_of_Includes>`_
      Use standard order for readability and to avoid hidden dependencies: Related header,
      C library, C++ library, other libraries' ``.h``, your project's ``.h``.

* **MC01##** Scoping

* **MC02##** Classes

    * **MC0202** (GCSG) `Copyable and Movable Types <https://google.github.io/styleguide/cppguide.html#Copyable_Movable_Types>`_
      Support copying and/or moving if it makes sense for your type. Otherwise, disable
      the implicitly generated special functions that perform copies and moves.

    * **MC0210** (GCSG) `Declaration Order <https://google.github.io/styleguide/cppguide.html#Declaration_Order>`_
      Use the specified order of declarations within a class: public: before private:,
      methods before data members (variables), etc.

* **MC03##** Functions

    * **MC0300** (GCSG) `Parameter Ordering <https://google.github.io/styleguide/cppguide.html#Function_Parameter_Ordering>`_

    * **MC0302** (GCSG) `Reference Arguments <https://google.github.io/styleguide/cppguide.html#Reference_Arguments>`_
      All parameters passed by reference must be labeled const.

* **MC04##** Other

    * **MC0401** (GCSG) `Variable-Length Arrays and alloca() <https://google.github.io/styleguide/cppguide.html#Variable-Length_Arrays_and_alloca__>`_
      We do not allow variable-length arrays or alloca().

    * **MC0407** (GCSG) `Preincrement and Predecrement <https://google.github.io/styleguide/cppguide.html#Preincrement_and_Predecrement>`_
      Use prefix form (``++i``) of the increment and decrement operators with iterators
      and other template objects.

    * **MC0408** (GCSG) `Use of const <https://google.github.io/styleguide/cppguide.html#Use_of_const>`_
      Use ``const`` whenever it makes sense. With C++11, ``constexpr`` is a better choice
      for some uses of ``const``.

* **MC05##** Naming

    * **MC0503** (GCSG) `Variable Names <https://google.github.io/styleguide/cppguide.html#Variable_Names>`_
      The names of variables and data members are all lowercase, with underscores between
      words. Data members of classes (but not structs) additionally have trailing
      underscores. For instance: ``a_local_variable``, ``a_struct_data_member``,
      ``a_class_data_member_``.

* **MC06##** Comments

* **MC07##** Formatting



API documentation
~~~~~~~~~~~~~~~~~

QA scripts will test if C++ source code documentation is missing.


Code coverage
~~~~~~~~~~~~~

For the moment, the C++ code is not included in the coverage analysis. It wasn't possible
to get ``gcov`` to work on the C++ extensions.


Mindless Cython criteria
========================

All Python criteria must be followed but nothing can be tested automatically. You have to
check everything manually. Because this is horribly inconvenient, the amount of Cython
code should be kept to a minimum.

In addition to the Python criteria, also use the following conventions:

* **MY0001** Pointers to NumPy array data should be accessed as follows:

  .. code-block:: python

    cdef double* pointer
    pointer = &array[0]         # for a 1D array
    pointer = &array[0, 0]      # for a 2D array
    pointer = &array[0, 0, 0]   # for a 3D array
    # etc.

* **MY0002** NumPy should be imported in Cython in a specific way. (See
  https://github.com/cython/cython/wiki/tutorials-numpy#c-api-initalization)

  .. code-block:: python

    import numpy as np
    cimport numpy as np
    np.import_array()
    # Other import and cimport lines should be put below.
