Git configuration and hooks
###########################

We recommend that you use the following ``~/.gitconfig`` file:

.. code-block:: ini

    [user]
        name = {Replace by your official name: First. M. Last}
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

Copy this script into the directory ``.git/hooks/pre-commit`` and make it
executable. The last part of the pre-commit script checks for python ``print``
lines. These should not be used in the Horton library. If you think you have
legitimate reasons to ignore this test, use the ``--no-verify`` option when
comitting.
