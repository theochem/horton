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

How to use Git
##############

HORTON uses `Git <http://git-scm.com/>`_ for version control.

A version control systems (VCS) allows many people to copy and modify the same source code while keeping
track of all changes made. The VCS software also helps you merge different
developments into one common source tree.

To refresh your mind on commonly used Git commands, please refer to `Git Reference <http://gitref.org/>`_.

This section goes through the basics of VCSs to get you started with developing new features in HORTON. Although the commands
below are specific to Git, the following entails good practices that can be
generalized and applied to other modern VCS such as Bazaar or Mercurial.


Installing Git
==============

The installation of git is covered in the sections :ref:`linux_install_dev` (Linux) or
:ref:`mac_install_dev` (Mac).


Git configuration
=================

We recommend you to set the following in your ``~/.gitconfig`` file:

.. code-block:: ini

    [user]
        name = {Replace by your official name: First M. Last}
        email = {Replace by a decent e-mail address, cfr. corresponding author on a paper.}

    [color]
        diff = auto
        status = auto
        interactive = auto
        branch = auto

    [push]
        default = simple

Also, install our pre-commit script as follows:

.. code-block:: bash

    cp -a tools/pre-commit .git/hooks/

This hook imposes some baseline quality checks on each commit:

.. literalinclude :: ../tools/pre-commit
    :language: bash
    :caption: tools/pre-commit

The last part of the ``pre-commit`` script checks for python ``print`` lines. These should
not be used in the HORTON library. If you think you have legitimate reasons to ignore this
check, use the ``--no-verify`` option when comitting.

Furthermore, it is useful to include the current branch in your shell prompt. To
do so, put one of the following in your ``~/.bashrc`` (Linux) or
``~/.bash_profile`` (Mac OS X) file:

* For terminals with a dark background:

  .. code-block:: bash

      GIT_PS='$(__git_ps1 ":%s")'
      export PS1="\[\033[1;32m\]\u@\h\[\033[00m\] \[\033[1;34m\]\w\[\033[00m\]\[\033[1;33m\]${GIT_PS}\[\033[1;34m\]>\[\033[00m\] "

* For terminals with a light background:

  .. code-block:: bash

      GIT_PS='$(__git_ps1 ":%s")'
      export PS1="\[\033[2;32m\]\u@\h\[\033[00m\]:\[\033[2;34m\]\w\[\033[3;31m\]${GIT_PS}\[\033[00m\]$ "

You can customize it to your taste. You may also want to add the ``export
PROMPT_DIRTRIM=3`` line to keep the shell prompt short. If you are a happy ``vim``
user, you can set ``export EDITOR=vim`` to get syntax highlighting when writing
commit messages.


Some terminology
================

Patch
    A set of changes in the source code. These are typically recorded in a
    `patch` file. Such a file specifies a set of lines that are removed and
    a set of lines that are added.

`SHA-1 <http://en.wikipedia.org/wiki/SHA-1>`_ hash
    A `numerical` checksum of a given length in bytes (in this case 256) for a
    much larger amount of data, e.g. a very long character string. There are usually
    two main goals when designing hashing algorithms: (i) it is not possible to
    derive the original data from a hash and (ii) a small change in the original
    data completely changes the hash. The `MD5
    <http://en.wikipedia.org/wiki/MD5>`_ checksum is well known and often used
    for CD images, but it is not great in terms of the above two hashing
    objectives.

Commit
    A patch with some extra information: author, timestamp, a SHA-1 hash of the
    code to which it applies, and some other things.

Branch
    A series of commits that describe the history of the source code.

    In realistic projects, the source code history is not linear, but contains
    many deviations from the ``master`` branch where people try to implement a
    new feature. It is, however, useful to have only one official linear history.
    We will show below how this can be done with git.

Branch HEAD
    The last commit in a branch.


Cloning the HORTON git repository
=================================

In order to `clone` the public HORTON repository, run the following commands:

.. code-block:: bash

    mkdir ~/code
    cd ~/code
    git clone git://github.com/theochem/horton.git
    cd horton

The version history can be updated with the latest committed patches on GitHub by:

.. code-block:: bash

    git pull

There is also a web interface to HORTON's git repository:
https://github.com/theochem/horton


.. _ref_build_refatoms:

Additional steps required to build the development version of HORTON
====================================================================

Several parts of HORTON make use of reference atomic computations. These files
are too large to be included in the git revision system. Therefore, they must be
downloaded separately when compiling a development version of HORTON:

.. code-block:: bash

    (cd data/refatoms; make all)


Work flow for adding a new feature
==================================

The development of a new feature typically consists of the following steps:

1. You make modifications of the code in a topic branch. You test and document your
   modifications, fix problems where needed.
2. Make a pull request on Github. (Some tests will be automatically executed.) Someone
   will review your pull request, which usually leads to suggestions to improve
   your modifications.
3. As soon as you pull request is up to snuff, it will be merged into the master branch.

.. note::

    Try to keep the amount of work in one branch as low as possible and get it
    reviewed/merged as early as possible. This takes some planning, as you have to
    figure out how to break your big plans up into smaller steps. In general
    this is a good exercise that will help you write more modular code.
    Although this seems to be cumbersome, it does save time for everyone involved.

When you intend to make relatively large modifications, it is recommended to discuss these
first, e.g. on the `HORTON mailing list
<https://groups.google.com/forum/#!forum/horton-discuss>`_, just to avoid disapointments
in the long run.


Develop the feature in a topic branch
---------------------------------------

0. `Fork <https://help.github.com/articles/fork-a-repo>`_ the public HORTON repository on
   Github (if not done yet), clone it on your local machine and enter the source tree:

   .. code-block:: bash

       $ ~/code> git clone https://github.com/your_account/horton.git
       $ ~/code> cd horton
       $ ~/.../horton:master>

   where ``your_account`` needs to be replaced by your Github account name.

1. Switch to the ``master`` branch, if needed:

   .. code-block:: bash

      $ ~/.../horton:foo> git checkout master
      $ ~/.../horton:master>

   Make sure there are no uncommitted changes in the source code on the ``foo``
   branch before switching to the ``master`` branch.

2. Get the latest version of the source code:

   .. code-block:: bash

    $ ~/.../horton:master> git pull origin

3. Make a topic branch, say ``bar``, and switch to it:

   .. code-block:: bash

    $ ~/.../horton:master> git checkout -b bar
    $ ~/.../horton:bar>

   Make sure that you are on the right branch before starting to implement the
   new feature ``bar``. (Try to pick a more meaningful branch name based on the feature
   you are implementing.)

4. Now you are in the right place to start making changes to the source code,
   and committing patches. When adding a new feature, also add
   tests, documentation, docstrings, comments and examples to clarify and debug the new feature.
   (The more tests, documentation and examples, the better.)

5. Review your changes with ``git diff``. Make sure there are no trailing white spaces
   or trailing blank lines. These can be removed with the ``./cleancode.sh``
   script. If you created new files, run the ``./updateheaders.py`` script to
   make sure the new files have the proper headers.

6. Get an overall overview of the added changes and new files with ``git status``.

7. Add the changed files that will be committed with ``git add <file_name>`` command. There are
   two ways to do this:

   * Add all changes in certain files:

     .. code-block:: bash

        $ ~/.../horton:bar> git add horton/file1.py horton/file2.py ...

   * Add interactively by going through the changes in all/some files:

     .. code-block:: bash

        $ ~/.../horton:bar> git add -p [horton/file1.py horton/file2.py ...]

8. Commit the added files to your working branch:

   .. code-block:: bash

      $ ~/.../horton:bar> git commit

   This command will start an editor in which you can write a commit message. By
   convention, such a message starts with a short single-line description
   of at most 69 characters. Optionally, a longer description follows
   that is separated from the short description by an empty line. More
   suggestions for writing meaningful commit messages can be found `here
   <http://chris.beams.io/posts/git-commit/>`_. If you only intend to write a
   short description, it can be included on the command line:

   .. code-block:: bash

      $ ~/.../horton:bar> git commit -m 'Short description'


In practice, you'll make a couple of commits before a new feature is finished. After
committing the changes and testing them thoroughly, you are ready for the next step.


Make your branch available for review with a pull request (PR)
--------------------------------------------------------------

In order to let others look at your code, you have to make your branch
available by pushing it to your forked Github repository.

1. Push your branch to the remote server:

   .. code-block:: bash

      git push origin bar:bar

2. Now go to the Github website and make a `Pull Request
   <https://help.github.com/articles/using-pull-requests/>`_ with the ``master`` branch of
   the ``theochem/horton`` repository as the destination. As soon as you do this, a series
   of basic QA tests will be executed to check for common problems. If these basic QA
   tests pass, someone will review your branch manually based on the
   :ref:`tech_dev_checklist`. You fix all the issues brought up during the review by
   making additional commits or, if you really messed up, by rewriting your branch. As
   soon as you push your changes back to the branch in your forked repository, they will
   show up in the PR, which triggers again the QA tests. When there are no further
   comments, your branch is ready to be merged.


Merging your pull request with the master branch
------------------------------------------------

You don't have to do anything for this, unless other branches got merged into
the master branch after you started your topic branch. In that case, you need to rebase
your topic branch on the current ``master`` branch and rerun all tests. This can be done
with the following steps:

1. `Synchronize <https://help.github.com/articles/syncing-a-fork/>`_ the ``master`` branch
   in your fork with the official HORTON repository.

3. Switch to your topic branch:

   .. code-block:: bash

      $ ~/.../horton:master> git checkout bar
      $ ~/.../horton:bar>

4. Create a new branch in which the result of ``git rebase`` will be stored:

   .. code-block:: bash

      $ ~/.../horton:bar> git checkout -b bar-1
      $ ~/.../horton:bar-1>


5. ``Rebase`` your commits on top of the latest ``master`` branch:

   .. code-block:: bash

      $ ~/.../horton:bar-1> git rebase master

   This command will try to apply the patches from your topic branch on top of the
   ``master`` branch. It may happen that changes in the ``master`` branch are not
   compatible with yours, such that your patches cannot be simply applied. When that is
   the case, the ``git rebase`` script will be interrupted and you are instructed on what
   to do. Do not panic when this happens. If you feel uncertain about how to resolve
   conflicts, it is time to call your git-savvy friends for help.

6. After the rebase procedure is complete, run all the tests again. If needed, fix
   problems and commit the changes.

7. Upload the commits to your fork:

   .. code-block:: bash

      $ ~/.../horton:bar-1> git push origin -f bar-1:bar

   This will rewrite the history of your topic branch, which will also show up in the PR.
   All automatic QA tests will be executed again.


Common issues
=============

* Remember to set the ``pre-commit`` hook. If this causes error messages when
  committing, use the ``cleancode.sh`` script. This removes all sorts of
  trailing white-space and converts every tab to four spaces. These conventions
  make ``git diff`` more meaningful and make it easier to merge and rebase commits.

* When you are customizing your bash prompt, you may get an error like
  ``__git_ps1: command not found...``, if you sourced ``git-completion.bash``.
  Then, before setting the ``GIT_PS``, you need to add the following line to your
  ``~/.bashrc`` (Linux) or ``~/.bash_profile`` (Mac OS X):

  .. code-block:: bash

     source /usr/share/git-core/contrib/completion/git-prompt.sh

  If you cannot find this file, you can get it from the link below:
  ``https://github.com/git/git/blob/master/contrib/completion/git-prompt.sh``
