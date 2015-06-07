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

How to use Git
##############

HORTON uses `Git <http://git-scm.com/>`_ for version control.

This section goes through the basics of version control. Although the commands
below are specific to Git, the following entails good practices that can be
generalized and applied to other modern version control systems (VCS) such
Bazaar or Mercurial.

A VCS allows many people to store and modify the same program while keeping
track of all the changes. VCS software also helps you merge different
developments into one common source tree.


Installing Git
==============

On most Linux distributions, git can be installed with a package manager:

* **Ubuntu Linux**:

  .. code-block:: bash

     sudo apt-get install git

* **Fedora Linux**:

  .. code-block:: bash

     sudo yum install git

On OS/X, an up-to-date version of git can be installed through MacPorts:

.. code-block:: bash

    sudo port install git


Git configuration
=================

We recommend that you use the following ``~/.gitconfig`` file:

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

The following ``pre-commit`` hook imposes some baseline quality checks on each
commit:

.. code-block:: bash

    #!/bin/bash

    red="\033[1;31m"
    color_end="\033[0m"

    # Check unwanted trailing whitespace or space/tab indents;
    if [[ `git diff --cached --check` ]]; then
        echo -e "${red}Commit failed: trailing whitespace, trailing empty lines, dos line endings${color_end}"
        git diff --cached --check
        exit 1
    fi

    # Check for untracked files (not in .gitignore)
    if [[ `git status -u data horton doc scripts tools -s | grep "^??"` ]]; then
        echo -e "${red}Commit failed: untracked files (not in .gitignore).${color_end}"
        git status -u data horton doc scripts tools -s | grep "^??"
        exit 1
    fi

    # Check for new print statements
    if [[ `git diff --cached | grep '^+' | sed  's/^.//' | sed 's:#.*$::g' | grep 'print '` ]]; then
        echo -e "${red}Commit failed: print statements${color_end}"
        git diff --cached | grep '^+' | sed  's/^.//' | sed 's:#.*$::g' | grep print
        exit 1
    fi

Copy this script into the directory ``.git/hooks/`` as ``pre-commit`` and make it
executable. The last part of the pre-commit script checks for python ``print``
lines. These should not be used in the HORTON library. If you think you have
legitimate reasons to ignore this test, use the ``--no-verify`` option when
comitting.

Furthermore, it is useful to include the current branch in your shell prompt. To
do so, put one of the following in your ``~/.bashrc`` (Linux) or
``~/.bash_profile`` (Max OS/X) file:

* For terminals with a dark background:

   .. code-block:: bash

      GIT_PS='$(__git_ps1 ":%s")'
      export PS1="\[\033[1;32m\]\u@\h\[\033[00m\] \[\033[1;34m\]\w\[\033[00m\]\[\033[1;33m\]${GIT_PS}\[\033[1;34m\]>\[\033[00m\] "

* For terminals with a light background:

   .. code-block:: bash

      GIT_PS='$(__git_ps1 ":%s")'
      export PS1="\[\033[2;32m\]\u@\h\[\033[00m\]:\[\033[2;34m\]\w\[\033[3;31m\]${GIT_PS}\[\033[00m\]$ "

You can customize it to your taste. You may also want to add a line ``export
PROMPT_DIRTRIM=3`` to keep the shell prompt short. If you are a happy ``vim``
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
    A patch with some extra information: author, timestamp, a SHA1 hash of the
    code to which it applies, and some other things.

Branch
    A series of commits that describe the history of the source code.

    In realistic projects, the source code history is not linear, but contains
    many deviations from the `master branch` where people try to implement a
    new feature. It is, however, useful to have only one official linear history.
    We will show below how this can be done with git.

Branch head
    The last commit in a branch.


Cloning the HORTON git repository
=================================

In order to `clone` the public HORTON repository, run these commands:

.. code-block:: bash

    mkdir ~/code
    cd ~/code
    git clone git://github.com/theochem/horton.git
    cd horton

The version history can be updated with the latest patches with the following
command:

.. code-block:: bash

    git pull

There is also a web interface to HORTON's git repository:
https://github.com/theochem/horton


Additional steps required to build the development version of HORTON
====================================================================

Several parts of HORTON make use of reference atomic computations. These files
are too large to be included in the git revision system. Therefore, they must be
downloaded separately when compiling a development version of HORTON:

.. code-block:: bash

    (cd data/refatoms; make all)


Work flow for adding a new feature
==================================

The development of a new feature typically consists of three large steps: (i)
modifications of the code in a separate branch, (ii) review of the new code,
fixing problems and (iii) rebase your branch on top of the `master` branch and
publish.

.. note::

    Try to keep the amount of work in one branch as low as possible and get it
    reviewed/merged as early as possible. This takes some planning, as you have to
    figure out how to break your big plans up into smaller steps. In general
    this is a good exercise that will help you to write more modular code.
    Although this seems to be a cumbersome approach, it does save time for
    everyone involved.

The instructions below are written for the general public, i.e. people that do
not have access the Clifford server. When you work with Clifford, the internal
development server, make the following substitutions below:

* ``master`` branch => ``prerelease`` branch
* Read-only ``origin`` repository at Github with URL ``https://github.com/theochem/horton.git`` =>
  Read-only ``origin`` repository with URL ``ssh://clifford/horton-release``
* Writable repository ``review`` for uploading your branches with a URL you created =>
  Writable repository ``review`` on clifford with URL ``ssh://clifford/horton-2``


Develop the feature in a separate branch
----------------------------------------

0. Clone the public HORTON repository (if not done yet) and enter the source
   tree:

   .. code-block:: bash

       $ ~/code> git clone https://github.com/theochem/horton.git
       $ ~/code> cd horton
       $ ~/.../horton:master>

1. Switch to the master branch if needed:

   .. code-block:: bash

      $ ~/.../horton:foo> git checkout master
      $ ~/.../horton:master>

   Make sure there are no uncommitted changes in the source code before
   switching to the master branch.

2. Get the latest version of the official code:

   .. code-block:: bash

    $ ~/.../horton:master> git pull origin

3. Make a new branch, e.g. named ``bar``:

   .. code-block:: bash

    $ ~/.../horton:master> git checkout -b bar
    $ ~/.../horton:bar>

   Only start changing the code and committing patches once you have changed
   to this new branch for the implementation of feature `bar`. (Try to pick
   a more meaningful branch name.)

4. Make some changes in the source code. When adding a new feature, also add
   tests, documentation, docstrings, comments and examples for that feature.
   (The more tests, documentation and examples, the better.)

5. Review your changes with ``git diff``. Make sure there are no trailing spaces
   or trailing blank lines. These can be removed with the ``./cleancode.sh``
   script. If you created new files, run the ``./updateheaders.py`` script to
   make sure the new files have the proper headers.

6. Review the changed/new files with ``git status``

7. Select the files/changes that will be committed with ``git add``. There are
   two ways to do this:

   * Add all changes in certain files:

     .. code-block:: bash

        $ ~/.../horton:bar> git add horton/file1.py horton/file2.py ...

   * Interactively go through the changes in all/some files:

     .. code-block:: bash

        $ ~/.../horton:bar> git add -p [horton/file1.py horton/file2.py ...]

8. Commit the selected files to your working branch:

   .. code-block:: bash

      $ ~/.../horton:bar> git commit

   This command will start an editor in which you can write a commit message. By
   convention, such a message starts with a short single-line description
   of at most 69 characters. Optionally, a longer description follows
   that is separated from the short description by an empty line. More
   suggestions for writing heplful commit messages can be found `here
   <http://chris.beams.io/posts/git-commit/>`_. If you only intend to write a
   short description, it can be provided at the command line:

   .. code-block:: bash

      $ ~/.../horton:bar> git commit -m 'Short description'


In practice, you'll make a few commits before a new feature is finished. After
adding a few commits, testing them thoroughly, you are ready for the next step.


Make your branch available for review
-------------------------------------

In order to let someone look at your code, you have to make your branch
available by pushing it to a remote server. One may use `Github
<http://www.github.com>`_ for this purpose.

1. Configure your repository for the remote server:

   .. code-block:: bash

      git remote add review <paste_your_remote_url_here>

2. Push your branch to the remote server:

   .. code-block:: bash

      git push review bar:bar

Now send the URL of your remote server and the name of the branch to a peer for
review. If you are looking for someone to review your code, post a request on
the `HORTON mailing list <https://groups.google.com/d/forum/horton-discuss>`_

Unless, you have written spotless code, you will make some further modifications
to the code, commit these and push them to the remote server for review. Once
this iterative process has converged, it is time to move to the next step.


Rebase your branch on top of the master branch
----------------------------------------------

It is likely that while developing your branch, the master branch
has evolved with new commits added by other developers. You need to append your
branch to the new HEAD of the master branch with ``git rebase``

1. Switch to the master branch:

   .. code-block:: bash

      $ ~/.../horton:bar> git checkout master
      $ ~/.../horton:master>

2. Get the latest version of the official code:

   .. code-block:: bash

      $ ~/.../horton:master> git fetch
      $ ~/.../horton:master> git pull

3. Switch to your working branch:

   .. code-block:: bash

      $ ~/.../horton:master> git checkout bar
      $ ~/.../horton:bar>

4. Create a new branch in which the result of ``git rebase`` will be stored.

   .. code-block:: bash

      $ ~/.../horton:bar> git checkout -b bar-1
      $ ~/.../horton:bar-1>


5. `Rebase` your commits on top of the latest master branch:

   .. code-block:: bash

      $ ~/.../horton:bar-1> git rebase master

    This command will try to apply the patches from your working branch to the
    master branch. It may happen that changes in the master branch are not
    compatible with yours, such that your patches cannot be simply applied.
    When that is the case, the ``git rebase`` script will interrupt and tell you
    what to do. Do not panic when this happens. If you feel uncertain about how
    to resolve conflicts, it is time to call your git-savvy friends for help.

6. After the rebase procedure is completed, run all tests again. If needed, fix
   problems and commit the changes.

7. Upload the commits to your remote server:

   .. code-block:: bash

      $ ~/.../horton:bar-1> git push review bar-1:bar-1

Now, you can get in touch with one of the HORTON developers (at the `HORTON
mailing list <https://groups.google.com/d/forum/horton-discuss>`_) to transfer
these new patches to the public master branch of HORTON.


Common Issues
=============

* Make sure you set the ``pre-commit`` hook. If this causes error messages when
  committing, use the script ``cleancode.sh``. This removes all sorts of
  trailing white-space and convert every tab to four spaces. These conventions
  make diffstats more meaningful and make it easier to merge and rebase commits.

* When you're customizing your bash prompt, you may get a bash error like
  ``__git_ps1: command not found...`` if you sourced ``git-completion.bash.bash``.
  Then, you need to add ``source /usr/share/git-core/contrib/completion/git-prompt.sh``
  before setting the ``GIT_PS``. If you cannot find this file, you can
  get it from ``https://github.com/git/git/blob/master/contrib/completion/git-prompt.sh``.
