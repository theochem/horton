#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Horton is a Density Functional Theory program.
# Copyright (C) 2011-2012 Toon Verstraelen <Toon.Verstraelen@UGent.be>
#
# This file is part of Horton.
#
# Horton is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
#
# Horton is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, see <http://www.gnu.org/licenses/>
#
#--


import os
import numpy as np
from glob import glob
from distutils.core import setup
from distutils.extension import Extension
from distutils.command.install_data import install_data
from Cython.Distutils import build_ext


def get_sources(dirname):
    # avoid accidental inclusion of in-place build files and inc files
    result = [fn for fn in glob('%s/*.cpp' % dirname)
              if not (('ext.cpp' in fn) or ('_inc.cpp' in fn))]
    result.append('%s/cext.pyx' % dirname)
    return result


def get_depends(dirname):
    result = glob('%s/*.h' % dirname)
    result += glob('%s/*.pxd' % dirname)
    result += glob('%s/*_inc.pxd' % dirname)
    return result




class my_install_data(install_data):
    """Add a datadir.txt file that points to the root for the data files. It is
       otherwise impossible to figure out the location of these data files at
       runtime.
    """
    def run(self):
        # Do the normal install_data
        install_data.run(self)
        # Create the file datadir.txt. It's exact content is only known
        # at installation time.
        dist = self.distribution
        libdir = dist.command_obj["install_lib"].install_dir
        for name in dist.packages:
            if '.' not in name:
                destination = os.path.join(libdir, name, "data_dir.txt")
                print "Creating %s" % destination
                if not self.dry_run:
                    f = file(destination, "w")
                    print >> f, self.install_dir
                    f.close()


setup(
    name='Horton',
    version='0.0',
    description='Horton is a Density Functional Theory program.',
    author='Toon Verstraelen',
    author_email='Toon.Verstraelen@UGent.be',
    url='http://molmod.ugent.be/code/',
    package_dir = {'horton': 'horton'},
    packages=['horton', 'horton.test',
              'horton.dpart', 'horton.dpart.test',
              'horton.grid', 'horton.grid.test',
              'horton.hamiltonian', 'horton.hamiltonian.test',
              'horton.gbasis', 'horton.gbasis.test'],
    cmdclass = {'build_ext': build_ext, 'install_data': my_install_data},
    data_files=[
        ('share/horton/', glob('data//*.*')),
        ('share/horton/test', glob('data/test/*.*')),
        ('share/horton/basis', glob('data/basis/*.*')),
        ('share/horton/refatoms', glob('data/refatoms/*.h5')),
    ],
    ext_modules=[
        Extension("horton.cext",
            sources=get_sources('horton'),
            depends=get_depends('horton'),
            include_dirs=[np.get_include()],
            language="c++"),
        Extension("horton.gbasis.cext",
            sources=get_sources('horton/gbasis'),
            depends=get_depends('horton/gbasis'),
            extra_objects=['depends/libint-2.0.0-stable/lib/libint2.a'],
            include_dirs=['depends/libint-2.0.0-stable/include', np.get_include()],
            language="c++"),
        Extension("horton.grid.cext",
            sources=get_sources('horton/grid'),
            depends=get_depends('horton/grid'),
            include_dirs=[np.get_include()],
            language="c++"),
        Extension("horton.hamiltonian.cext",
            sources=get_sources('horton/hamiltonian'),
            depends=['depends/libxc-1.2.0/src/.libs/libxc.a'],
            extra_objects=['depends/libxc-1.2.0/src/.libs/libxc.a'],
            include_dirs=['depends/libxc-1.2.0/src', np.get_include()],
            language="c++"),
    ],
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License (GPL)',
        'Operating System :: POSIX :: Linux',
        'Programming Language :: Python',
        'Topic :: Science/Engineering :: Molecular Science'
    ],
)
