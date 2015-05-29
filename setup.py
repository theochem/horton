#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Horton is a development platform for electronic structure methods.
# Copyright (C) 2011-2015 The Horton Development Team
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


import os, sys, platform, ConfigParser, subprocess
import numpy as np
from glob import glob
from distutils.core import setup
from distutils.extension import Extension
from distutils.command.install_data import install_data
from distutils.command.install_headers import install_headers
from Cython.Distutils import build_ext


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
        return subprocess.check_output(['pkg-config', libname, '--'+option],
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
    return all(len(value)==0 for value in lib_config.itervalues())


def all_exist(lib_config):
    '''Test if all paths in the lib_config exist'''
    for key, value in lib_config.iteritems():
        for path in value:
            if not os.path.exists(path):
                return False
    return True


def detect_machine():
    '''Return a description of the machine name, used for data/setup_cfgs/...'''
    if sys.platform == 'linux2':
        dist = platform.linux_distribution()
        return ('Linux-%s-%s-%s' % (dist[0], dist[1], platform.machine())).replace(' ', '_')
    elif sys.platform == 'darwin':
        mac_ver = platform.mac_ver()
        return 'Darwin-%s-%s' % (mac_ver[0], mac_ver[2])
    else:
        return 'unknown'


def lib_config_magic(prefix, libname, static_config={}):
    '''Detect the configuration of a given library

       **Arguments:**

       prefix
            The prefix for this library. This is a name that Horton uses to
            refer to the library.

       libname
            The library name as it is known to the compiler and to pkg-config.
            For example, if the shared object is libfoo.so, then the library
            name is foo.

       **Optional arguments**

       static_config
            If given, this static library configuration is attempted. Ignored
            when empty, or when it contains non-existing files.
    '''
    print '%s Configuration' % prefix.upper()

    # Start out empty
    lib_config = dict((key, []) for key in lib_config_keys)

    # Update with info from setup.cfg
    lib_config.update(get_lib_config_setup(prefix, 'setup.cfg'))

    # Override with environment variables
    lib_config.update(get_lib_config_env(prefix))

    # If no environment variables were set, attempt to guess the right settings.
    if all_empty(lib_config):
        if all_exist(static_config) and not all_empty(static_config):
            # If the static build is present, use it.
            lib_config.update(static_config)
        else:
            try:
                # Try to get dynamic link info from pkg-config
                lib_config.update(get_lib_config_pkg(libname))
            except PkgConfigError:
                # Try to find a machine-specific config file
                machine = detect_machine()
                fn_setup_cfg = 'data/setup_cfgs/setup.%s.cfg' % machine
                lib_config.update(get_lib_config_setup(prefix, fn_setup_cfg))
                if all_empty(lib_config):
                    # Uber-dumb fall back
                    lib_config['libraries'] = [libname]

    print_lib_config('Final', lib_config)
    return lib_config


# Print the Machine name on screen
# --------------------------------

print 'MACHINE=%s' % detect_machine()


# Configuration of LibXC
# ----------------------

# Static build info in the depends directory to check for:
libxc_dir = 'depends/libxc-2.2.2'
libxc_static_config = {
    'extra_objects': ['%s/src/.libs/libxc.a' % libxc_dir],
    'include_dirs': ['%s/src' % libxc_dir, libxc_dir],
}
# Detect the configuration for LibXC
libxc_config = lib_config_magic('libxc', 'xc', libxc_static_config)


# Configuration of LibInt2
# ------------------------

libint2_dir = 'depends/libint-2.0.3-stable'
libint2_static_config = {
    'extra_objects': ['%s/lib/.libs/libint2.a' % libint2_dir],
    'include_dirs': ['%s/include' % libint2_dir],
}
libint2_config = lib_config_magic('libint2', 'int2', libint2_static_config)


# Configuration of BLAS
# ---------------------

# First try to get the BLAS Configuration from environment variables.
blas_config = lib_config_magic('blas', 'atlas')

# Detect which BLAS implementation is used and set a corresponding preprocessor
# flag.
blas_names = blas_config['libraries'] + blas_config['extra_objects']
if any(('mkl' in l) for l in blas_names):
    blas_precompiler = ('BLAS_MKL', '1')
elif any(('atlas' in l) for l in blas_names):
    blas_precompiler = ('BLAS_ATLAS', '1')
elif any(('openblas' in l) for l in blas_names):
    blas_precompiler = ('BLAS_OPENBLAS', '1')
else:
    print '   Unknown BLAS implementation. Assuming Netlib-compatible headers.'
    blas_precompiler = ('BLAS_OTHER', '1')
print 'BLAS precompiler directive: -D%s' % blas_precompiler[0]


# Call distutils setup
# --------------------

setup(
    name='horton',
    version='2.0.0',
    description='Horton is a development platform for electronic structure methods.',
    author='Toon Verstraelen',
    author_email='Toon.Verstraelen@UGent.be',
    url='http://theochem.github.com/horton/',
    scripts=glob("scripts/*.py"),
    package_dir = {'horton': 'horton'},
    packages=['horton', 'horton.test',
              'horton.correlatedwfn', 'horton.correlatedwfn.test',
              'horton.espfit', 'horton.espfit.test',
              'horton.gbasis', 'horton.gbasis.test',
              'horton.grid', 'horton.grid.test',
              'horton.io', 'horton.io.test',
              'horton.localization', #'horton.localization.test',
              'horton.matrix', 'horton.matrix.test',
              'horton.meanfield', 'horton.meanfield.test',
              'horton.orbitalentanglement', #'horton.orbitalentanglement.test',
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
            include_dirs=[np.get_include(), '.'] +
                          libint2_config['include_dirs'] +
                          blas_config['include_dirs'],
            library_dirs=libint2_config['library_dirs'] +
                         blas_config['library_dirs'],
            libraries=libint2_config['libraries'] + blas_config['libraries'],
            extra_objects=libint2_config['extra_objects'] +
                          blas_config['extra_objects'],
            extra_compile_args=libint2_config['extra_compile_args'] +
                                blas_config['extra_compile_args'],
            extra_link_args=libint2_config['extra_link_args'] +
                             blas_config['extra_link_args'],
            define_macros=[blas_precompiler],
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
            include_dirs=[np.get_include(), '.'] + libxc_config['include_dirs'],
            library_dirs=libxc_config['library_dirs'],
            libraries=libxc_config['libraries'],
            extra_objects=libxc_config['extra_objects'],
            extra_compile_args=libxc_config['extra_compile_args'],
            extra_link_args=libxc_config['extra_link_args'],
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
