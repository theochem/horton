#!/usr/bin/env python
# -*- coding: utf-8 -*-
# HORTON: Helpful Open-source Research TOol for N-fermion systems.
# Copyright (C) 2011-2017 The HORTON Development Team
#
# This file is part of HORTON.
#
# HORTON is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
#
# HORTON is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, see <http://www.gnu.org/licenses/>
#
# --


import ConfigParser
import distutils.ccompiler
from distutils.command.install_data import install_data
from distutils.command.install_headers import install_headers
from distutils.core import setup
from distutils.extension import Extension
from glob import glob
import json
import os
import platform
import subprocess
import sys

import numpy as np
from Cython.Distutils import build_ext


# Distutils optimizations
# -----------------------

def parallelCCompile(self, sources, output_dir=None, macros=None,
                     include_dirs=None, debug=0, extra_preargs=None, extra_postargs=None,
                     depends=None):
    """Monkey-patch for distutils compiler to run in parallel."""
    # these lines are copied from distutils.ccompiler.CCompiler directly
    macros, objects, extra_postargs, pp_opts, build = self._setup_compile(
        output_dir, macros, include_dirs, sources, depends, extra_postargs)
    cc_args = self._get_cc_args(pp_opts, debug, extra_preargs)
    # parallel code
    N = 2  # number of parallel compilations
    import multiprocessing.pool
    def _single_compile(obj):
        try:
            src, ext = build[obj]
        except KeyError:
            return
        self._compile(obj, src, ext, cc_args, extra_postargs, pp_opts)

    # convert to list, imap is evaluated on-demand
    list(multiprocessing.pool.ThreadPool(N).imap(_single_compile, objects))
    return objects

distutils.ccompiler.CCompiler.compile = parallelCCompile


# Utility functions
# -----------------

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
        # at installation time. By default, it is the installation prefix
        # passed to setup.py, but one can override it using the env var
        # INSTALL_DIR, which may be useful for packaging, or any other
        # situation where the installed files are moved to a new location
        # afterwards.
        my_install_dir = os.getenv("INSTALL_DIR", self.install_dir)
        # Loop over all packages in this project and write the data_dir.txt
        # file only in the main package. Usualy, there is only one that matters.
        dist = self.distribution
        libdir = dist.command_obj["install_lib"].install_dir
        for name in dist.packages:
            # If a package contains a dot, e.g. horton.test, then don't write
            # the file data_dir.txt.
            if '.' not in name:
                destination = os.path.join(libdir, name, "data_dir.txt")
                print "install_dir={}".format(my_install_dir)
                print "Creating {}".format(destination)
                if not self.dry_run:
                    with open(destination, "w") as f:
                        print >> f, my_install_dir


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


# Library configuration functions
# -------------------------------

lib_config_keys = ['include_dirs', 'library_dirs', 'libraries', 'extra_objects',
                   'extra_compile_args', 'extra_link_args']


def print_lib_config(heading, lib_config):
    '''Print (partial) lib_config'''
    print '   %s' % heading
    if len(lib_config) == 0:
        print '      -'
    else:
        for key, value in sorted(lib_config.iteritems()):
            if len(value) > 0:
                print '      %s: %s' % (key, value)


def get_lib_config_setup(prefix, fn_setup_cfg):
    '''Get library configuration from a setup.cfg'''
    lib_config = {}
    if os.path.isfile(fn_setup_cfg):
        config = ConfigParser.ConfigParser()
        config.read(fn_setup_cfg)
        if config.has_section(prefix):
            for key in lib_config_keys:
                if config.has_option(prefix, key):
                    value = config.get(prefix, key).strip()
                    if value is not None and len(value) > 0:
                        lib_config[key] = value.split(':')
        print_lib_config('From %s' % fn_setup_cfg, lib_config)
    else:
        print '   File %s not found. Skipping.' % fn_setup_cfg
    return lib_config


def get_lib_config_env(prefix):
    '''Read library config from the environment variables'''
    lib_config = {}
    for key in lib_config_keys:
        varname = ('%s_%s' % (prefix, key)).upper()
        value = os.getenv(varname)
        if value is not None:
            lib_config[key] = value.split(':')
    print_lib_config('From environment variables', lib_config)
    return lib_config


class PkgConfigError(Exception):
    pass


def run_pkg_config(libname, option):
    '''Safely try to call pkg-config'''
    try:
        return subprocess.check_output(['pkg-config', libname, '--' + option],
                                       stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError:
        raise PkgConfigError('pkg-config did not exit properly')
    except OSError:
        raise PkgConfigError('pkg-config not installed')


def get_lib_config_pkg(libname):
    '''Get library config from the pkg-config program'''
    lib_config = {
        'include_dirs': [word[2:] for word in run_pkg_config(libname, 'cflags-only-I').split()],
        'library_dirs': [word[2:] for word in run_pkg_config(libname, 'libs-only-L').split()],
        'libraries': [word[2:] for word in run_pkg_config(libname, 'libs-only-l').split()],
        'extra_compile_args': run_pkg_config(libname, 'cflags-only-other').split(),
        'extra_link_args': run_pkg_config(libname, 'libs-only-other').split(),
    }
    print_lib_config('From pkg-config', lib_config)
    return lib_config


def all_empty(lib_config):
    '''Test if all lib_config fields are empty'''
    if len(lib_config) == 0:
        return True
    return all(len(value) == 0 for value in lib_config.itervalues())


def all_exist(lib_config):
    '''Test if all paths in the lib_config exist'''
    for key, value in lib_config.iteritems():
        for path in value:
            if not os.path.exists(path):
                return False
    return True


def lib_config_magic(prefix, libname, static_config={}, known_include_dirs=[]):
    '''Detect the configuration of a given library

    Parameters
    ----------

    prefix : str
        The prefix for this library. This is a name that HORTON uses to refer to the
        library.

    libname : str
        The library name as it is known to the compiler and to pkg-config. For example, if
        the shared object is libfoo.so, then the library name is foo.

    static_config : dict
        If given, this static library configuration is attempted. Ignored when empty, or
        when it contains non-existing files.

    known_include_dirs : list of str
        When all other methods of finding the library settings fail, the first existing
        directory in this list is added to the include path. This is useful when header
        files are commonly installed in a place that is not considered by default by most
        compilers.
    '''
    print '%s Configuration' % prefix.upper()

    # Start out empty
    lib_config = dict((key, []) for key in lib_config_keys)

    # Update with info from setup.cfg
    lib_config.update(get_lib_config_setup(prefix, 'setup.cfg'))

    # Override with environment variables
    lib_config.update(get_lib_config_env(prefix))

    # If no environment variables were set, attempt to use the static config.
    if all_empty(lib_config):
        if all_empty(static_config):
            print '   No static config available for this library'
        elif not all_exist(static_config):
            print_lib_config('Static lib not found in ${QAWORKDIR}', static_config)
        else:
            # If the static build is present, use it.
            print_lib_config('Static lib config in ${QAWORKDIR}', static_config)
            lib_config.update(static_config)

    # If also the static config did not work, try pkg-config
    if all_empty(lib_config):
        try:
            # Try to get dynamic link info from pkg-config
            lib_config.update(get_lib_config_pkg(libname))
        except PkgConfigError:
            print '   pkg-config failed.'

    # Uber-dumb fallback. It works most of the times.
    if all_empty(lib_config):
        lib_config['libraries'] = [libname]
        for include_dir in known_include_dirs:
            if os.path.isdir(include_dir):
                lib_config['include_dirs'] = [include_dir]
                break
        print_lib_config('Last resort fallback plan', lib_config)

    print_lib_config('Final', lib_config)
    return lib_config


# Print the Machine name on screen
# --------------------------------

print 'PLATFORM={}'.format(platform.platform())

# Load dependency information
# ---------------------------
with open('dependencies.json') as f:
    dependencies = json.load(f)
# Order does not matter here. Just make it easy to look things up
dependencies = dict((d['name'], d) for d in dependencies)

# Locate ${QAWORKDIR}
# -------------------
qaworkdir = os.getenv('QAWORKDIR')
if qaworkdir is None:
    qaworkdir = 'qaworkdir'

# Configuration of LibXC
# ----------------------

# Static build info in the QAWORKDIR:
libxc_dir = '%s/cached/libxc-%s' % (qaworkdir, str(dependencies['libxc']['version_ci']))
libxc_static_config = {
    'extra_objects': ['%s/lib/libxc.a' % libxc_dir],
    'include_dirs': ['%s/include' % libxc_dir],
}
# Common include dirs that are not considered by the compiler by default:
known_libxc_include_dirs = ['/opt/local/include']
# Detect the configuration for LibXC
libxc_config = lib_config_magic(
    'libxc', 'xc', libxc_static_config, known_libxc_include_dirs)

# Configuration of LibInt2
# ------------------------

# Static build info in the QAWORKDIR:
libint2_dir = '%s/cached/libint-%s' % (qaworkdir, str(dependencies['libint']['version_ci']))
libint2_static_config = {
    'extra_objects': ['%s/lib/libint2.a' % libint2_dir],
    'include_dirs': ['%s/include/libint2' % libint2_dir],
}
# Common include dirs that are not considered by the compiler by default:
known_libint2_include_dirs = ['/usr/include/libint2', '/opt/local/include/libint2']
libint2_config = lib_config_magic(
    'libint2', 'int2', libint2_static_config, known_libint2_include_dirs)

# Print versions of (almost) all dependencies
# -------------------------------------------
print 'Version of dependencies:'
for name, dependency in sorted(dependencies.iteritems()):
    version_command = dependency.get('version_command')
    if version_command is not None:
        try:
            version_info = subprocess.check_output(
                dependency['version_command'], shell=True,
                stderr=subprocess.STDOUT).strip()
        except subprocess.CalledProcessError:
            version_info = '-- not found --'
        print '{:>20}: {}'.format(name, version_info)


# Define extension modules
# ------------------------

ext_modules = [
    Extension(
        "horton.cext",
        sources=get_sources('horton'),
        depends=get_depends('horton'),
        include_dirs=[np.get_include(), '.'],
        extra_compile_args=['-std=c++11'],
        language="c++"),
    Extension(
        "horton.gbasis.cext",
        sources=get_sources('horton/gbasis'),
        depends=get_depends('horton/gbasis'),
        include_dirs=[np.get_include(), '.'] +
                     libint2_config['include_dirs'],
        library_dirs=libint2_config['library_dirs'],
        libraries=libint2_config['libraries'],
        extra_objects=libint2_config['extra_objects'],
        extra_compile_args=libint2_config['extra_compile_args'] +
                           ['-std=c++11'],
        extra_link_args=libint2_config['extra_link_args'],
        language="c++"),
    Extension(
        "horton.grid.cext",
        sources=get_sources('horton/grid') + [
              'horton/cell.cpp',
              'horton/moments.cpp'],
        depends=get_depends('horton/grid') + [
              'horton/cell.pxd', 'horton/cell.h',
              'horton/moments.pxd', 'horton/moments.h'],
        include_dirs=[np.get_include(), '.'],
        extra_compile_args=['-std=c++11'],
        language="c++", ),
    Extension(
        "horton.meanfield.cext",
        sources=get_sources('horton/meanfield'),
        depends=get_depends('horton/meanfield'),
        include_dirs=[np.get_include(), '.'] + libxc_config['include_dirs'],
        library_dirs=libxc_config['library_dirs'],
        libraries=libxc_config['libraries'],
        extra_objects=libxc_config['extra_objects'],
        extra_compile_args=libxc_config['extra_compile_args'] + ['-std=c++11'],
        extra_link_args=libxc_config['extra_link_args'],
        language="c++"),
    Extension(
        "horton.espfit.cext",
        sources=get_sources('horton/espfit') + [
            'horton/cell.cpp',
            'horton/grid/uniform.cpp'],
        depends=get_depends('horton/espfit') + [
            'horton/cell.pxd', 'horton/cell.h',
            'horton/grid/uniform.pxd', 'horton/grid/uniform.h'],
        include_dirs=[np.get_include(), '.'],
        extra_compile_args=['-std=c++11'],
        language="c++"),
]

for e in ext_modules:
    e.cython_directives = {"embedsignature": True}

# Call distutils setup
# --------------------

setup(
    name='horton',
    version='2.1.0',
    description='HORTON: Helpful Open-source Research TOol for N-fermion systems.',
    author='Toon Verstraelen',
    author_email='Toon.Verstraelen@UGent.be',
    url='http://theochem.github.com/horton/',
    scripts=glob("scripts/*.py"),
    package_dir={'horton': 'horton'},
    packages=[
        'horton', 'horton.test',
        'horton.espfit', 'horton.espfit.test',
        'horton.gbasis', 'horton.gbasis.test',
        'horton.grid', 'horton.grid.test',
        'horton.io', 'horton.io.test',
        'horton.meanfield', 'horton.meanfield.test',
        'horton.part', 'horton.part.test',
        'horton.scripts', 'horton.scripts.test',
        'horton.modelhamiltonians', 'horton.modelhamiltonians.test'],
    cmdclass={
        'build_ext': build_ext,
        'install_data': my_install_data,
        'install_headers': my_install_headers,
    },
    data_files=[
        ('share/horton', glob('data/*.*')),
        ('share/horton/test', glob('data/test/*.*')),
        ('share/horton/basis', glob('data/basis/*.*')),
        ('share/horton/grids', glob('data/grids/*.txt')),
        ('share/horton/refatoms', glob('data/refatoms/*.h5')),
    ] + [
        ('share/horton/examples/%s' % os.path.basename(dn[:-1]),
         glob('%s/*.py' % dn) + glob('%s/README' % dn))
        for dn in glob('data/examples/*/')
    ] + [
        ('include/horton', glob('horton/*.h')),
        ('include/horton/grid', glob('horton/grid/*.h')),
        ('include/horton/gbasis', glob('horton/gbasis/*.h')),
        ('include/horton/espfit', glob('horton/espfit/*.h')),
    ],
    package_data={
        'horton': ['*.pxd'],
        'horton.espfit': ['*.pxd'],
        'horton.gbasis': ['*.pxd'],
        'horton.grid': ['*.pxd'],
    },
    ext_modules=ext_modules,
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
    ], requires=['numpy', 'nose']
)
