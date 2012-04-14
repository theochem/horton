Getting Started with Coding
###########################

This chapter discusses all the technology and conventions that make it fun to
work on Horton with a small team, and to make sure things remain manageable
as the elephant grows bigger and bigger.


Coding guidelines
=================

Although the development of Horton is a rather informal business, it is not a
bad idea to follow some guidelines to guarantee the quality of the program. It
may feel a bit frustrating in the beginning, but it avoids a lot of frustration
in the long run. These guidelines assume you are familiar with the revision
system, with Python in general, and know how to test to code. (vide infra)

1. **Write validation routines for new pieces of code. Commit new code and tests
   in one patch.**

   Writing tests is difficult, but really necessary. Only fools believe they can
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
   It is also OK to upload these branches to the central repository, such that
   others can see and review your code. As soon as one of your branches is
   mature enough, use ``git rebase`` to apply your patches to the master branch.

4. **When you find a bug, first write a test that fails because of the bug. Then
   fix the bug. Commit both test and fix in one patch.**

   This is a great way of fixing bugs and making tests. It also guarantees that
   the bug you found will not be reintroduced accidentally in later versions of
   Horton.

5. **Write lots of documentation, doc-strings and comments.**

   `Comments` are needed for obvious reasons.

   `Doc-strings` are a rather Pythonic thing to do. More information about
   doc-strings can be found `here <http://www.python.org/dev/peps/pep-0257#what-is-a-docstring>`_.

   `The documentation` is written in `ReStructuredText <http://docutils.sourceforge.net/rst.html>`_
   and compiled with `Sphinx <http://sphinx.pocoo.org/>`_ into the fancy
   webpage you are currently looking at. More details on the documentation
   system are given below.

6. **Adopt a clean coding style.**

   There is an extensive `Style Guide for Python Code <http://www.python.org/dev/peps/pep-0008/>`_.
   The main purpose is to make the source code easier to read and maintain.
   These are the basics:

   * Never use tabs. Always use four spaces to indent code.
   * Use CamelCase for class names. All other names are written in lower case,
     with underscores to separate different words in a name.
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


Git: the fast version control system
====================================

Horton uses `Git <http://git-scm.com/>`_ for version control.

This section goes through
the basics. Although the commands below are specific for Git, most of the
good practices are general and also apply to other modern VCSs such Bazaar
or Mercurial. There are however a few old-fashioned systems such as CVS (and to
some extent Subversion), that do not (or only poorly) support some of the
guidelines below. Of all the VCSs out there, Git is probably the most snappy, and
definitely the most popular for new projects. Git is also the only widely
adopted open-source VCS for the development of commercial software.

The purpose of a VCS is to allow many different people to make modifications in the same
program without messing things up. VCS software also makes life much easier for
someone who has to merge different developments into one common source tree.
With the help of a few good practices, it becomes even unnecessary to have a
codemeister at all.


Configuring Git
---------------

Make a file ~/.gitconfig as follows::

    [user]
        name = Toon Verstraelen
        email = Toon.Verstraelen@UGent.be

    [color]
        diff = always
        status = always
        interactive = always
        branch = always

    [branch]
        autosetuprebase = always

    [branch "master"]
        remote = origin
        merge = master

    [push]
        default = matching

Replace my name and email by yours. If you are also working on other projects
that use git, it may be useful to move some of these options to the file
``.git/config`` in the Horton source tree. Add the following to your
``.bashrc``, and modify to your taste::

    GITPS1='$(__git_ps1 ":%s")'
    SILVER="\[\033[0;37m\]"
    GRAY="\[\033[1;30m\]"
    GREEN="\[\033[1;32m\]"
    BLUE="\[\033[1;34m\]"
    YELLOW="\[\033[1;33m\]"
    RS="\[\033[00m\]"

    export PROMPT_DIRTRIM=3
    export PS1="${GREEN}\u@\h${RS} ${BLUE}\w${RS}${YELLOW}${GITPS1}${BLUE}>${RS} "

With these settings, the bash prompt will show what `branch` you are working on.
It will become clear later what that means.


Some terminology
----------------

Patch
    A set of changes in the source code. These are typically recorded in a
    `patch` file. Such a file specifies a set of lines that are removed and
    a set of lines that are added.

`SHA-1 <http://en.wikipedia.org/wiki/SHA-1>`_ hash
    A `numerical` checksum of a given length in bytes (in this case 256) for a
    much larger amount of data, e.g. a very long character string. One tries to
    design hashing algorithms such that they are doing two things very well: (i) it
    is not possible to derive the original data from a hash and (ii) a small
    change in the original data completely changes the hash. The `MD5 <http://en.wikipedia.org/wiki/MD5>`_ checksum is
    well known and often used from CD images, but it is not great in terms of
    the above two hashing objectives.

Commit (git specific)
    A patch with a some extra information: author, timestamp, a SHA1 hash of the
    code to which it applies, and some other things.

Branch
    A series of commits that describe the history of the source code.

    In realistic projects, the source code history is not linear, but contains
    many deviations from the `official branch` where people try to implement a
    new feature. It is however useful to have only one official linear history.
    We will show below how this can be done with git.

Branch head
    The last commit in a branch.


Adding a feature
----------------

Only the basic work flow is discussed, so things may become more complicated.

1. Switch to the master branch if needed::

    toony@poony ~/.../horton:foo> git checkout master
    toony@poony ~/.../horton:master>

   The master branch is the official branch of Horton. Also make sure there
   are no uncommitted changes in the source code before switching to the
   master branch.

2. Get the latest version of the official code::

    toony@poony ~/.../horton:master> git pull

3. Make a new branch::

    toony@poony ~/.../horton:master> git checkout -b bar
    toony@poony ~/.../horton:bar>

   Only start changing the code and committing patches once you have changed
   to a dedicated branch for the implementation of feature `bar`.

4. Make some changes in the source code. When adding a new feature, also add
   tests for that feature. (The more tests, the better.)

5. Review your changes with ``git diff``. Make sure there are no trailing spaces
   or trailing blank lines. These can be removed with the ``./cleancode.sh``
   script. If you created new files, run the ``./updateheaders.py`` script to
   make sure the new files have the proper headers.

6. Review the changed/new files with ``git status``

7. Select the files/changes that will be committed with ``git add``. There are
   two ways to do this:

   * Add all changes in certain files::

        toony@poony ~/.../horton:bar> git add horton/file1.py horton/file2.py ...

   * Interactively go through the changes in all/some files::

        toony@poony ~/.../horton:bar> git add -p [horton/file1.py horton/file2.py ...]

8. Commit the selected files to your working branch::

    toony@poony ~/.../horton:bar> git commit -m 'Short description'

In practice, you'll make a few commits before a new feature is finished. After
adding a few commits, testing them thoroughly and having them reviewed by a
peer, you may want to transfer your commits from your working branch ``bar`` to
the master branch. This is done as follows:

1. Switch to the master branch::

    toony@poony ~/.../horton:bar> git checkout master
    toony@poony ~/.../horton:master>

2. Get the latest version of the official code::

    toony@poony ~/.../horton:master> git pull

3. Switch to your working branch::

    toony@poony ~/.../horton:master> git checkout bar
    toony@poony ~/.../horton:bar>

4. `Rebase` your commits on top of the latest master branch::

    toony@poony ~/.../horton:bar> git rebase master

..

    This command will try to apply the patches from your working branch to the
    master branch. It may happen that others have changed the official version
    such that your patches do no longer simply apply. When that is the case,
    the ``git rebase`` script will interrupt and tell you what to do.

5. Run all tests again once the rebase procedure is completed.

6. Upload the commits::

    toony@poony ~/.../horton:bar> git push origin bar:master

   This may fail when someone has uploaded yet a few more patches in the
   meantime. Start again at step 1 to update your patches.

The last step will only work if you have write access to the main Horton
repository. If this is not the case, there are a few other ways to get your
contributions into the master branch. The best option is to upload your patches
to a personal git repository and notify one of the authors with write access
about your contributions. They should take care of reviewing and committing your
patches.


Git shortcuts
-------------

The above procedures are rather lengthly, and for small changes there are
shortcuts that are more convenient.

1. One may commit all changed files (not new ones) without using ``git add`` first::

    git commit -a -m 'Short description'

2. One does not have to switch to a working branch to commit a few patches. This
   can also be done in the master branch with less commands:

   a. Get the latest official patches (in the master branch)::

        toony@poony ~/.../horton:master> git pull

   b. Make some changes, test and commit::

        toony@poony ~/.../horton:master> git commit -a -m 'Short description'

   c. Do a combined rebase and pull::

        toony@poony ~/.../horton:master> git pull --rebase

      If someone has added patches to the master branch in the meantime, they
      will be downloaded. Your patches are then applied afterwards with ``git
      rebase`` automatically.

   d. Upload your commits::

        toony@poony ~/.../horton:master> git push

   The downside of this approach is that one may accidentally make merges if
   one does not carefully follow the instructions. Merges make the history
   of the master branch non-linear, which is messy and hard to follow.


Writing tests
=============

Horton uses the `Nosetests <http://somethingaboutorange.com/mrl/projects/nose/0.11.2/>`_
program to run all validation routines. Use one of the existing tests as an
example, or go through the Nosetests documentation to learn how to write tests
from scratch.

All tests in Horton are located in the directory ``horton/test``. All module
files containing tests have a filename that starts with ``test_``. In these
modules, all functions with a name that starts with ``test_`` are picked up
by Nosetests. Tests that do not follow this convention, are simply ignored.

The tests are run as follows (including preparation steps)::

    toony@poony ~/.../horton:master> ./cleanfiles.sh
    toony@poony ~/.../horton:master> ./setup.py build_ext -i
    toony@poony ~/.../horton:master> nosetests -v

There are some cases where the first two commands are not needed. You will
figure out.


Writing documentation
=====================

All the documentation is located in the ``doc`` directory. We use the `Sphinx
<http://sphinx.pocoo.org/>`_ formatting engine to compile the `documentation
source code` into fancy formatted HTML or PDF.

The source files have the extension ``.rst``, and are written in the
`ReStructuredText <http://docutils.sourceforge.net/rst.html>`_ (RST) format.
RST in some sense comparable to latex, but more intuitive to use.
It also has some specific advantages for documenting software.

All ``.rst``-files are part of the source tree, just like the actual source
code. Git is also used to keep track of changes in the documentation.

There is a makefile to generate the documentation based in the source code::

    toony@poony ~/.../horton:master> cd doc
    toony@poony ~/.../horton/doc:master> make html
    toony@poony ~/.../horton/doc:master> make pdf


Add yourself as an author
=========================

The authors of Horton are listed in a few places. This matters for the
open-source license. (In the long run, it would be better to have a small
consortium that can hold the copyrights, but it is too early for that.)

* If you are going to write source code, update the source headers and commit a
  patch as follows:

  1. Edit HEADER and add your name and email
  2. Run the script ``./updateheaders.py``
  3. ``git commit -a -m 'Hello horton source code!'``

* If you are going to write documentation, update the file ``doc/conf.py`` and
  commit that change. There is no need to run some update script.
