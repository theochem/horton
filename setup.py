#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Horton is a development platform for electronic structure methods.
# Copyright (C) 2011-2013 Toon Verstraelen <Toon.Verstraelen@UGent.be>
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



libintdir = 'depends/libint-2.0.3-stable'
libxcdir = 'depends/libxc-2.0.3'
execfile('customize.py')


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



# Configure the linkage of the extension that use libint and libxc. If the
# libararies are compiled in the depends directory, these are used for
# static linking. The user may also change the directories in which the
# hand-compiled libraries are found, i.e. by making the necessary changes
# in customize.py Otherwhise, dynamic linking is attempted.

if os.path.isfile('%s/src/.libs/libxc.a' % libxcdir):
    libxc_extra_objects = ['%s/src/.libs/libxc.a' % libxcdir]
    libxc_include_dirs = ['%s/src' % libxcdir, libxcdir]
    libcx_libraries = []
else:
    libxc_extra_objects = []
    libxc_include_dirs = []
    libcx_libraries = ['xc']

if os.path.isfile('%s/lib/.libs/libint2.a' % libintdir):
    libint_extra_objects = ['%s/lib/.libs/libint2.a' % libintdir]
    libint_include_dirs = ['%s/include' % libintdir]
    libint_libraries = []
else:
    libint_extra_objects = []
    libint_include_dirs = []
    libint_libraries = ['int2']


setup(
    name='horton',
    version='1.2.1',
    description='Horton is a development platform for electronic structure methods.',
    author='Toon Verstraelen',
    author_email='Toon.Verstraelen@UGent.be',
    url='http://theochem.github.com/horton/',
    scripts=glob("scripts/*.py"),
    package_dir = {'horton': 'horton'},
    packages=['horton', 'horton.test',
              'horton.espfit', 'horton.espfit.test',
              'horton.gbasis', 'horton.gbasis.test',
              'horton.grid', 'horton.grid.test',
              'horton.io', 'horton.io.test',
              'horton.meanfield', 'horton.meanfield.test',
              'horton.part', 'horton.part.test',
              'horton.scripts', 'horton.scripts.test',
              'horton.modelhamiltonians',
              'horton.modelhamiltonians.test'],
    cmdclass = {'build_ext': build_ext, 'install_data': my_install_data},
    data_files=[
        ('share/horton/', glob('data/*.*')),
        ('share/horton/test', glob('data/test/*.*')),
        ('share/horton/basis', glob('data/basis/*.*')),
        ('share/horton/grids', glob('data/grids/*.txt')),
        ('share/horton/refatoms', glob('data/refatoms/*.h5')),
    ] + [
        ('share/horton/examples/%s' % os.path.basename(dn[:-1]), glob('%s/*.py' % dn) + glob('%s/*.xyz' % dn))
        for dn in glob('data/examples/*/')
    ],
    ext_modules=[
        Extension("horton.cext",
            sources=get_sources('horton'),
            depends=get_depends('horton'),
            include_dirs=[np.get_include()],
            language="c++"),
        Extension("horton.gbasis.cext",
            sources=get_sources('horton/gbasis') + ['horton/moments.cpp'],
            depends=get_depends('horton/gbasis') + ['horton/moments.pxd', 'horton/moments.h'],
            extra_objects=libint_extra_objects,
            libraries=libint_libraries,
            include_dirs=[np.get_include(), 'horton'] + libint_include_dirs,
            language="c++"),
        Extension("horton.grid.cext",
            sources=get_sources('horton/grid') + [
                'horton/cell.cpp',
                'horton/moments.cpp'],
            depends=get_depends('horton/grid') + [
                'horton/cell.pxd', 'horton/cell.h',
                'horton/moments.pxd', 'horton/moments.h'],
            include_dirs=[np.get_include(), 'horton'],
            #extra_compile_args=["-fopenmp"],
            #extra_link_args=["-fopenmp"],
            language="c++",),
        Extension("horton.meanfield.cext",
            sources=get_sources('horton/meanfield'),
            depends=get_depends('horton/meanfield'),
            extra_objects=libxc_extra_objects,
            libraries=libcx_libraries,
            include_dirs=[np.get_include()] + libxc_include_dirs,
            language="c++"),
        Extension("horton.espfit.cext",
            sources=get_sources('horton/espfit') + [
                'horton/cell.cpp',
                'horton/grid/uniform.cpp'],
            depends=get_depends('horton/espfit') + [
                'horton/cell.pxd', 'horton/cell.h',
                'horton/grid/uniform.pxd', 'horton/grid/uniform.h'],
            include_dirs=[np.get_include(), 'horton', 'horton/grid'],
            language="c++"),
    ],
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License (GPL)',
        'Operating System :: POSIX :: Linux',
        'Programming Language :: Python :: 2',
        'Programming Language :: Cython',
        'Programming Language :: C++',
        'Topic :: Science/Engineering :: Molecular Science'
    ],
)
