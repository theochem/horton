Making Pull requests
####################

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
        default = simple

Replace my name and email by yours. If you are also working on other projects
that use git, it may be useful to move some of these options to the file
``.git/config`` in the Horton source tree.


Furthermore, it is useful to include the current branch in your shell promt. Put
one of the following in your ``~/.bashrc`` file:

* For terminals with a dark background::

    GIT_PS='$(__git_ps1 ":%s")'
    export PS1="\[\033[1;32m\]\u@\h\[\033[00m\] \[\033[1;34m\]\w\[\033[00m\]\[\033[1;33m\]${GIT_PS}\[\033[1;34m\]>\[\033[00m\] "

* For terminals with a light background::

    GIT_PS='$(__git_ps1 ":%s")'
    export PS1="\[\033[2;32m\]\u@\h\[\033[00m\]:\[\033[2;34m\]\w\[\033[3;31m\]${GIT_PS}\[\033[00m\]$ "

Add salt and pepper to taste. You may also want to add a line ``export
PROMPT_DIRTRIM=3`` to keep the shell prompt short.


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


Work flow for adding a new feature
----------------------------------

The development of a new feature typically consists of three large steps: (i)
modifications to the code in a separate branch, (ii) review of the new code,
fixing problems and (iii) rebase your branch on top of the master branch and
publish.


Develop the feature in a separate branch
........................................

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
   to this new branch for the implementation of feature `bar`.

4. Make some changes in the source code. When adding a new feature, also add
   tests, documentation, docstrings, comments and examples for that feature.
   (The more tests and documentation, the better.)

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
adding a few commits, testing them thoroughly, you are ready for the next step.


Make your branch available for review
.....................................

In order to let someone look at your code, you have to make your branch
available by pushing it to a remote server. One may use `Github
<http://www.github.com>`_ for this purpose.

1. Configure your repository for the remote server::

    git remote add review <paste_your_remote_url_here>

2. Push your branch to the remote server::

    git push remote bar:bar

Now send the URL of your remote server and the name of the branch to a peer for
review. If you are looking for someone to review your code, post a request on
the `Horton mailing list <https://groups.google.com/d/forum/horton-discuss>`_

Unless, you have written spotless code, you will make some further modifications
to the code, commit these and push them to the remote server for review. Once
this iterative process has converged, it is time to move to the next step.


Rebase your branch on top of the master branch
..............................................

It is likely that during the development of your feature, the master branch
has evolved with new commits added by other developers. You need to append your
branch to the new HEAD of the master branch with the program ``git rebase``

1. Switch to the master branch::

    toony@poony ~/.../horton:bar> git checkout master
    toony@poony ~/.../horton:master>

2. Get the latest version of the official code::

    toony@poony ~/.../horton:master> git pull

3. Switch to your working branch::

    toony@poony ~/.../horton:master> git checkout bar
    toony@poony ~/.../horton:bar>

4. Create a new branch in which the result of ``git rebase`` will be stored.

    toony@poony ~/.../horton:bar> git checkout -b bar-1
    toony@poony ~/.../horton:bar-1>


4. `Rebase` your commits on top of the latest master branch::

    toony@poony ~/.../horton:bar-1> git rebase master

..

    This command will try to apply the patches from your working branch to the
    master branch. It may happen that others have changed the official version
    such that your patches do no longer simply apply. When that is the case,
    the ``git rebase`` script will interrupt and tell you what to do. Do not
    panic when this happens. If you feel uncertain about how to resolve
    conflicts, it is time to call your git-savvy friends for help.

5. Run all tests again once the rebase procedure is completed. If needed fix
   problems and commit the changes.

6. Upload the commits to your remote server::

    toony@poony ~/.../horton:bar-1> git push review bar-1:bar-1

Now you can get in touch with one of the Horton developers (at the `Horton
mailing list <https://groups.google.com/d/forum/horton-discuss>`_) to transfer
these new patches to the public master branch of Horton.
