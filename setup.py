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
from distutils.command.install_headers import install_headers
from Cython.Distutils import build_ext



libintdir = 'depends/libint-2.0.3-stable'
libxcdir = 'depends/libxc-2.1.0'
execfile('customize.py')


def get_sources(dirname):
    '''Get all cpp files and the cext.pyx file of a package'''
    # avoid accidental inclusion of in-place build files and inc files
    result = [fn for fn in glob('%s/*.cpp' % dirname)
              if not (('ext.cpp' in fn) or ('_inc.cpp' in fn))]
    result.append('%s/cext.pyx' % dirname)
    return result


def get_depends(dirname):
    '''Get all files that should trigger a recompilation of the C extension of a package'''
    result = glob('%s/*.h' % dirname)
    result += glob('%s/*.pxd' % dirname)
    return result


def get_headers():
    '''Get all header-like files that need to be installed'''
    result = []
    for dn in ['horton/'] + glob('horton/*/'):
        result.extend(glob('%s/*.h' % dn))
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


class my_install_headers(install_headers):
    def run(self):
        headers = self.distribution.headers
        if not headers:
            return

        self.mkpath(self.install_dir)
        for header in headers:
            dest = os.path.join(os.path.dirname(self.install_dir), header)
            dest_dn = os.path.dirname(dest)
            if not os.path.isdir(dest_dn):
                self.mkpath(dest_dn)
            (out, _) = self.copy_file(header, dest)
            self.outfiles.append(out)



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
    libint_include_dirs = ['/usr/include/libint2']
    libint_libraries = ['int2']

# Switches to MKL BLAS from default ATLAS BLAS. Specify with environmental variable
# USE_MKL=1. You are responsible for setting the CPATH and configuring ldconfig yourself.
if os.getenv("USE_MKL") == "1":
    blas_include_dirs = []
    blas_lib_dirs = []
    blas_libraries = ['mkl_intel_lp64', 'mkl_core', 'mkl_sequential', 'pthread',
            'm']
    blas_precompiler = ("CBLAS_CPP", "<mkl.h>")
else:
    blas_include_dirs = ["/usr/include/atlas-x86_64-sse3"]
    blas_lib_dirs = ['/usr/lib64/atlas-sse3']
#    blas_lib_dirs = []
    blas_libraries = ['atlas','cblas']
    blas_precompiler = ("CBLAS_C", "<cblas.h>")



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
              'horton.correlatedwfn',
              'horton.espfit', 'horton.espfit.test',
              'horton.gbasis', 'horton.gbasis.test',
              'horton.grid', 'horton.grid.test',
              'horton.io', 'horton.io.test',
              'horton.matrix', 'horton.matrix.test',
              'horton.meanfield', 'horton.meanfield.test',
              'horton.part', 'horton.part.test',
              'horton.scripts', 'horton.scripts.test',
              'horton.modelhamiltonians',
              'horton.modelhamiltonians.test'],
    cmdclass = {
        'build_ext': build_ext,
        'install_data': my_install_data,
        'install_headers': my_install_headers,
    },
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
    package_data={
        'horton': ['*.pxd'],
        'horton.espfit': ['*.pxd'],
        'horton.gbasis': ['*.pxd'],
        'horton.grid': ['*.pxd'],
    },
    ext_modules=[
        Extension("horton.cext",
            sources=get_sources('horton'),
            depends=get_depends('horton'),
            include_dirs=[np.get_include(), '.'],
            language="c++"),
        Extension("horton.matrix.cext",
            sources=get_sources('horton/matrix'),
            depends=get_depends('horton/matrix'),
            include_dirs=[np.get_include(), '.'],
            language="c++"),
        Extension("horton.gbasis.cext",
            sources=get_sources('horton/gbasis') + ['horton/moments.cpp'],
            depends=get_depends('horton/gbasis') + ['horton/moments.pxd', 'horton/moments.h'],
            extra_objects=libint_extra_objects,
            libraries=libint_libraries+blas_libraries,
            include_dirs=[np.get_include(), '.'] + libint_include_dirs +
            blas_include_dirs,
            library_dirs = blas_lib_dirs,
            define_macros = [blas_precompiler],
            language="c++"),
        Extension("horton.grid.cext",
            sources=get_sources('horton/grid') + [
                'horton/cell.cpp',
                'horton/moments.cpp'],
            depends=get_depends('horton/grid') + [
                'horton/cell.pxd', 'horton/cell.h',
                'horton/moments.pxd', 'horton/moments.h'],
            include_dirs=[np.get_include(), '.'],
            #extra_compile_args=["-fopenmp"],
            #extra_link_args=["-fopenmp"],
            language="c++",),
        Extension("horton.meanfield.cext",
            sources=get_sources('horton/meanfield'),
            depends=get_depends('horton/meanfield'),
            extra_objects=libxc_extra_objects,
            libraries=libcx_libraries,
            include_dirs=[np.get_include(), '.'] + libxc_include_dirs,
            language="c++"),
        Extension("horton.espfit.cext",
            sources=get_sources('horton/espfit') + [
                'horton/cell.cpp',
                'horton/grid/uniform.cpp'],
            depends=get_depends('horton/espfit') + [
                'horton/cell.pxd', 'horton/cell.h',
                'horton/grid/uniform.pxd', 'horton/grid/uniform.h'],
            include_dirs=[np.get_include(), '.'],
            language="c++"),
    ],
    headers=get_headers(),
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
