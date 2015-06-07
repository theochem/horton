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

Common problems
###############

Trouble building or installing HORTON
=====================================

Let's assume you have already built and installed all your dependencies correctly.
You try to install HORTON, i.e.

.. code-block:: bash

    ./setup.py install --user

or you run nosetests, but something goes wrong. Remember, you have installed all your dependencies
correctly, so the problems you are facing has to do with finding and using those
dependencies. Chances are, you have installed the dependencies using different
methods, e.g. yum, apt-get, macports, homebrew, pip, and compilation from source. You
may even have the same dependency installed multiple times using different methods.
You have to make sure your system can find the right dependency and be able to use
it. We've encountered problems with three types of dependencies: Python modules,
executables, and libraries. Also, multiple installation of Python may cause problems.

Python modules
--------------

If you have installed some python module (i.e. numpy, scipy, cython, h5py, sympy,
matplotlib, nosetests, sphinx, doxygen, breathe, docutils) and you get an error
saying your system cannot find such module, then you need to set the environment
variable ``PYTHONPATH`` to include your module's location. By default, python
already includes some paths without having to set ``PYTHONPATH``. These directories are
stored in attribute ``path`` of module ``sys``, which can be accessed by:

.. code-block:: bash

    python -c "import sys; import pprint; pprint.pprint(sys.path)"

You should make sure that there are no modules of the same name in all of these
directories. If there are recurrences, then you should make sure that you are
selecting the module that you want by ensuring that the desired directory comes before
the others. Please note that this is quite a bad practice, and will likely
result in more confusion in the future. In the following examples, we may
use directories that are already included in Python by default. For some
special cases, these directories were excluded, and needed to be included manually.

For example, some Mac users needed to set
the ``PYTHONPATH`` after instaslling modules through PIP:

.. code-block:: bash

    export PYTHONPATH=${HOME}/Library/Python/2.7/lib/python/site-packages:$PYTHONPATH

or their system site-packages:

.. code-block:: bash

    export PYTHONPATH=/Library/Python/2.7/lib/python/site-packages:$PYTHONPATH

Similarly, a few Linux users needed to set ``PYTHONPATH`` after installation through PIP:

.. code-block:: bash

    export PYTHONPATH=${HOME}/.local/lib/python2.7/site-packages

or

.. code-block:: bash

    export PYTHONPATH=/lib/python2.7/site-packages

or

.. code-block:: bash

    export PYTHONPATH=/lib64/python2.7/site-packages


The idea is to find your module (the directory that contains the files to use that
module), and add its parent directory to the ``PYTHONPATH``. Finding the module
can be tricky, but remember, there are finite number of files on your hard drive.
You can use unix commands like ``find`` to speed up the process. You can look up
where your package manager installs packages. You can explicitly build your module
in a location that you will remember (try to adopt some convention for installing
modules). Including the python modules are commonly installed in

.. code-block:: bash

   ${HOME}/Library/Python/2.7/site-packages
   /Library/Python/2.7/site-packages

for Mac users, and

.. code-block:: bash

   ${HOME}/.local/lib/python2.7/site-packages
   /lib/python2.7/site-packages
   /lib64/python2.7/site-packages

for Linux users.


Excecutables
------------

Let's say for whatever reason, HORTON requires the use of an executable for
installation. Then, this executable must be in the same directory as the current
directory, i.e. root of your source tree, or in one of the directories in ``PATH``
environment variable. By default, bash already includes some paths without
having to set ``PATH``. The contents of ``PATH`` can be accessed by:

.. code-block:: bash

    echo $PATH

For some special cases, these directories were excluded, often because ``PATH``
is overwritten somewhere in the ``${HOME}/.profile`` or ``${HOME}/.bashrc`` for
Mac and linux users, respectively. You should make sure that there are no
executables of the same name in all of the directories in ``PATH``. If there
are recurrences, then you should make sure that you are selecting the executable
that you want by ensuring that the desired directory comes before the others.
Please note that this is quite a bad practice, and will likely result in more
confusion in the future. In the following examples, we may use directories
that are already included in ``PATH`` by default. For some cases these directories
need to be added again to ``PATH``.

For example, Mac users that uses python scripts might do

.. code-block:: bash

    export PATH=${HOME}/Library/Python/2.7/bin:$PATH
    export PATH=/Library/Python/2.7/bin:$PATH

Similarly, Linux users may do

.. code-block:: bash

    export PATH=${HOME}/.local/bin:$PATH
    export PATH=/bin:$PATH

Using linux function ``find`` may help you find the appropriate directory.

Libraries
---------

Let's assume you have built your library correctly. Then, you need to make sure
``setup.py`` can find your libraries and their executables. You should consult
:ref:`setup_cfg` for a more complete understanding of the library linking process
in HORTON installation. Here, we will show how we solved some problems we've faced
in finding and linking libraries.

First, we need to find the library. We can locate the libraries install by using
the unix command ``ldconfig``:

.. code-block:: bash

    ldconfig -p | grep libraryname

``ldconfig -p`` prints all cached libraries, and piping to ``grep`` searches through
the results for the library with the ``libraryname``. Perhaps ``find`` would give
a more thorough search, especially if your library has not been cached yet. Here is
an example, where we tried to find atlas libraries in a cluster:

.. code-block:: bash

    ldconfig -p | grep atlas

which gives

.. code-block:: bash

	libptf77blas.so.3 (libc6,x86-64) => /usr/lib64/atlas/libptf77blas.so.3
	libptf77blas.so (libc6,x86-64) => /usr/lib64/atlas/libptf77blas.so
	libptcblas.so.3 (libc6,x86-64) => /usr/lib64/atlas/libptcblas.so.3
	libptcblas.so (libc6,x86-64) => /usr/lib64/atlas/libptcblas.so
	liblapack.so.3 (libc6,x86-64) => /usr/lib64/atlas/liblapack.so.3
	liblapack.so (libc6,x86-64) => /usr/lib64/atlas/liblapack.so
	libf77blas.so.3 (libc6,x86-64) => /usr/lib64/atlas/libf77blas.so.3
	libf77blas.so (libc6,x86-64) => /usr/lib64/atlas/libf77blas.so
	libclapack.so.3 (libc6,x86-64) => /usr/lib64/atlas/libclapack.so.3
	libclapack.so (libc6,x86-64) => /usr/lib64/atlas/libclapack.so
	libcblas.so.3 (libc6,x86-64) => /usr/lib64/atlas/libcblas.so.3
	libcblas.so (libc6,x86-64) => /usr/lib64/atlas/libcblas.so
	libatlas.so.3 (libc6,x86-64) => /usr/lib64/atlas/libatlas.so.3
	libatlas.so (libc6,x86-64) => /usr/lib64/atlas/libatlas.so

We see that all the libraries are located in ``/usr/lib64/atlas/``. Notice that
all the libraries are in x86-64 instruction set.

Then, we need to find the include directory. You can find this using the ``find``
function. Usually, the include directory is almost same as the library directory,
except instead of the ``lib``, there is ``include``. Continuing the above example,

.. code-block:: bash

    ls -d /usr/include/*atlas*

will give the list of directories that includes the word ``atlas``. The output
gives:

.. code-block:: bash

   /usr/include/atlas
   /usr/include/atlas-x86_64-base

Since we used the x86-64 instruction set, we select the directory that would
correspond with that instruction set, i.e. ``/usr/include/atlas-x86_64-base``.

In the above list of libraries associated with atlas, we have ``ptf77blas``,
``ptcblas``, ``lapack``, ``f77blas``, ``clapack``, ``cblas``, and ``atlas``.
Though we can include all these libraries, HORTON only uses ``atlas`` and
``cblas``. Therefore, the resulting ``setup.cfg`` file includes

.. code-block:: bash

  [blas]
  library_dirs=/usr/lib64/atlas
  libraries=atlas:cblas
  include_dirs=/usr/include/atlas-x86_64-base


Similarly, we can repeat the process for the LibXC and Libint2, where the libraries
that are needed are only ``libxc`` and ``libint``, respectively. See :ref:`setup_cfg`
for more details.

Multiple Python
---------------

Ideally, you should only have one Python in your system. Having multiple versions
may result in conflicts when running Python. You can use the ``find`` function to
see where the other Python's may be installed. First, you should find the location of
Python you are currently using:

.. code-block:: bash

    which python

If this directory does not correspond to the Python you want to use, you may want to
rename the binary file, and create a symbolic link for the Python you want to use.
For example, if you have Python 3 installed, and it is set as the default Python,
you may want to set Python 2.7 as the default:

.. code-block:: bash

    export temp=which python
    cd ${temp}
    rm python
    ln -s python2.7 python
    rm python-config
    ln -s python2.7-config python-config

where files ``python`` and ``python2.7-config`` were removed because it was a soft
link to ``python3.3``. If it was not a soft link, you may want to save this file
by renaming it, i.e.

.. code-block:: bash

    mv python python3.3_2
    mv python-config python3.3-config_2

You can check if a file is a link by the unix command

.. code-block:: bash

    ls -l

which gives output such as

.. code-block:: bash

    lrwxrwxrwx. 1 root root     9 Feb 25 04:28 python3 -> python3.3
    -rwxr-xr-x. 2 root root 11296 Dec  4  2014 python3.3

where the first file was a link and the second was not.
