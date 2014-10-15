Coding guidelines
#################

The following guidelines may feel a bit frustrating in the beginning, but it
avoids a lot of frustration in the long run.

1. **Write unit tests (validation routines) for new pieces of code. Commit new
   code and tests in one patch.**

   Writing unit tests is difficult, especially for newcomers, but it is really
   necessary to maintain to keep Horton functional. Only fools believe they can
   write bug-free code without testing. It always takes some creativity to
   devise test routines, certainly in Horton. Some examples are given below, but
   feel free to blur the lines. Just make sure a test does not take too long to
   execute (i.e. not more than a second).

   * Check whether a computational result is translational invariant.
   * Compare energies from atomic computations with published results, or
     results obtained with other (reliable) programs.
   * Compare analytical derivatives with finite differences.
   * Compare the output of a routine with a results derived on paper.

   In most cases, it is sufficient to imagine what `could` go wrong if your new
   code `would` be buggy.

   *Important:* Avoid tests that just check the output of a routine with the output
   from a previous version of the code. Such tests need to be redone each time
   a bug is found in that routine.

2. **Never commit a patch that breaks the tests.**

   The purpose of the tests is to keep the program working as we add more
   features. When a patch breaks previous tests, it is simply unacceptable.

3. **Ask someone to review your patches before they go into the official master
   branch.**

   It is perfectly OK to work in your own branch and do there whatever you like.
   It is also OK to upload these branches, such that others can see and review
   your code. As soon as one of your branches is mature enough, use ``git
   rebase`` to apply your patches to the master branch.

4. **When you find a bug, first write a test that fails because of the bug. Then
   fix the bug. Commit both test and fix in one patch.**

   This is a great way of fixing bugs and making tests. It also guarantees that
   the bug you found will not be reintroduced accidentally in later versions of
   Horton.

5. **Write lots of comments, doc-strings, documentation and examples.**

   `Comments` are needed for obvious reasons.

   `Doc-strings` are a rather Pythonic thing to do. More information about
   doc-strings can be found `here <http://www.python.org/dev/peps/pep-0257#what-is-a-docstring>`_.

   `The documentation` is written in `ReStructuredText <http://docutils.sourceforge.net/rst.html>`_
   and compiled with `Sphinx <http://sphinx.pocoo.org/>`_ into the fancy
   webpage you are currently looking at. More details on the documentation
   system are given below.

   `Examples` are simple python scripts that use Horton to run a computation.
   Keep these examples simple, e.g. no more than 20 lines, and add abundant
   comments that explain every line.

6. **Adopt a clean coding style.**

   There is an extensive `Style Guide for Python Code <http://www.python.org/dev/peps/pep-0008/>`_.
   The main purpose is to make the source code easier to read and maintain.
   These are the basics:

   * Never use tabs. Always use four spaces to indent code.
   * Use UpperCamelCase for class names. All other names are written in lower
     case, with underscores to separate different words in a name. Never use
     headlessCamelCase.
   * Separate top-level blocks in a module by two blank lines.
   * Separate methods in a class with one blank line.
   * Doc-strings also have a specific format:
        1. A one-line summary, followed by a blank line.
        2. Description of the arguments, followed by a blank line.
        3. Description of the return value(s), followed by a blank line.
           (optional)
        4. Some extra paragraphs with detailed information. (optional)
   * Use self-explaining names for classes, variables, functions, methods. Start
     functions and methods with a verb.
   * Never use ``from foo import *`` in a Horton module. Always specify what is
     to be imported, or use the ``import foo`` or ``import foo as bar`` syntax.
   * Avoid long lines, i.e. no longer than 80 characters. Use the line-wrapping
     features in Python to break long lines into smaller ones.

   Specific conventions for Horton

   * Import the Numpy package as follows: ``import numpy as np``.
   * Import the H5Py package as follows: ``import h5py as h5``.
   * Use single-letter variable names ``i``, ``j``, ``k`` and ``l`` only for
     integer loop variables.
