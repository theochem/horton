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


import glob
from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

setup(
    name='Horton',
    version='0.0',
    description='Horton is a Density Functional Theory program.',
    author='Toon Verstraelen',
    author_email='Toon.Verstraelen@UGent.be',
    url='http://molmod.ugent.be/code/',
    package_dir = {'horton': 'horton'},
    packages=['horton', 'horton/test',
              'horton/grid', 'horton/grid/test',
              'horton/gbasis', 'horton/gbasis/test'],
    cmdclass = {'build_ext': build_ext},
    ext_modules=[
        Extension("horton.gbasis.cext",
            sources=['horton/gbasis/boys.cpp', 'horton/gbasis/cext.pyx',
                     'horton/gbasis/cartpure.cpp', 'horton/gbasis/common.cpp',
                     'horton/gbasis/gbasis.cpp', 'horton/gbasis/ints.cpp',
                     'horton/gbasis/iter_gb.cpp', 'horton/gbasis/iter_pow.cpp'],
            depends=['horton/gbasis/boys.h', 'horton/gbasis/boys_inc.cpp', 'horton/gbasis/boys.pxd',
                     'horton/gbasis/cartpure.h', 'horton/gbasis/cartpure.pxd',
                     'horton/gbasis/common.h', 'horton/gbasis/common.pxd'
                     'horton/gbasis/gbasis.h', 'horton/gbasis/gbasis.pxd',
                     'horton/gbasis/ints.h', 'horton/gbasis/ints.pxd',
                     'horton/gbasis/iter_gb.h', 'horton/gbasis/iter_gb.pxd',
                     'horton/gbasis/iter_pow.h', 'horton/gbasis/iter_pow.pxd'],
            extra_objects=['depends/libint-2.0.0-stable/lib/libint2.a'],
            include_dirs=['depends/libint-2.0.0-stable/include'],
            language="c++"),
        Extension("horton.grid.cext",
            sources=['horton/grid/becke.cpp', 'horton/grid/cext.pyx',
                     'horton/grid/lebedev_laikov.cpp'],
            depends=['horton/grid/becke.h', 'horton/grid/becke.pxd',
                     'horton/grid/lebedev_laikov.h', 'horton/grid/lebedev_laikov.pxd'],
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
