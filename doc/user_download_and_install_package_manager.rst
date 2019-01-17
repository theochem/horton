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


Download and installation with a package manager
################################################

The installation of HORTON is supported by some package managers, which may be a lot
easier than the from-source installation discussed in the following sections. Below we try
to keep track of them.

Installing HORTON with any of these methods requires less technical skills but it may have
disadvantages too: (i) less control over the compilation process and (ii) older versions
of HORTON and its dependencies.

After installation with one of the package managers, you should be able to run the tests
as follows:

.. code-block:: bash

    nosetests -v horton


MacPorts
========

MacPorts is a package manager specific for OSX. After you have installed MacPorts, see
https://www.macports.org/install.php, you can install HORTON as follows:

.. code-block:: bash

    sudo port install gcc49; sudo port select --set gcc mp-gcc49
    sudo port install python27 +readline; sudo port select --set python python27
    sudo port install horton
    sudo port install py27-pip; sudo port select --set pip pip27

This will also install all the dependencies required for HORTON.

At the time of writing, there is still a small glitch in the MacPorts installation of the
LibXC dependency, which may cause some tests to fail.

There seems to be no official way to use MacPorts without ``sudo``, while this may
potentially break your system. You've been warned.


Easybuild
=========

EasyBuild is a package manager for HPC environments but works on any Linux system. Once
you have EasyBuild installed and configured, see http://hpcugent.github.io/easybuild/, you
can install HORTON e.g. as follows:

.. code-block:: bash

   eb horton-2.0.0-intel-2015b-Python-2.7.10-HDF5-1.8.15-patch1-serial.eb --robot

The ``--robot`` option will instruct EasyBuild to install all dependencies as well.

There are several config files one can use for different version of HORTON, in combination
with different compilers and Python dependencies. For the current list of configs, see:
https://github.com/hpcugent/easybuild-easyconfigs/tree/master/easybuild/easyconfigs/h/horton
You can change the toolchain and the version of HORTON with the ``--try-*`` options
discussed here:
http://easybuild.readthedocs.io/en/latest/Using_the_EasyBuild_command_line.html#tweaking-existing-easyconfig-files-using-try
