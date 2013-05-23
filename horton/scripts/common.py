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



import os, sys, datetime, numpy as np, h5py as h5, time, contextlib

from horton import UniformGrid, angstrom, periodic, Cell, log, dump_hdf5_low


__all__ = [
    'iter_elements', 'reduce_data', 'parse_h5', 'parse_ewald_args', 'parse_pbc',
    'parse_ugrid', 'store_args', 'safe_open_h5', 'write_part_output',
]


def iter_elements(elements_str):
    '''Interpret a string as a list of elements

       elements_str
            A string with comma-separated element numbers. One may add ranges
            with the format 'N-M' where M>N.
    '''
    for item in elements_str.split(','):
        if '-' in item:
            words = item.split("-")
            if len(words) != 2:
                raise ValueError("Each item should contain at most one dash.")
            first = periodic[words[0]].number
            last = periodic[words[1]].number
            if first > last:
                raise ValueError('first=%i > last=%i' % (first, last))
            for number in xrange(first,last+1):
                yield number
        else:
            yield periodic[item].number


def reduce_data(cube_data, ugrid, stride, chop):
    '''Reduce the uniform grid data according to stride and chop arguments

       **Arguments:**

       fn_cube
            The cube file with the electron density.

       ugrid
            The uniform integration grid.

       stride
            The reduction factor.

       chop
            The number of slices to chop of the grid in each direction.

       Returns: a new array and an updated grid object
    '''
    if (chop < 0):
        raise ValueError('Chop must be positive or zero.')
    if ((ugrid.shape - chop) % stride != 0).any():
        raise ValueError('The stride is not commensurate with all three grid demsions.')


    if chop == 0:
        new_cube_data = cube_data[::stride, ::stride, ::stride].copy()
    else:
        new_cube_data = cube_data[:-chop:stride, :-chop:stride, :-chop:stride].copy()

    new_shape = (ugrid.shape-chop)/stride
    grid_rvecs = ugrid.grid_cell.rvecs*stride
    new_ugrid = UniformGrid(ugrid.origin, grid_rvecs, new_shape, ugrid.pbc)

    return new_cube_data, new_ugrid


def parse_h5(arg_h5):
    if arg_h5.count(':') != 1:
        raise ValueError('An hdf5 argument must contain one colon.')
    return arg_h5.split(':')


def parse_ewald_args(args):
    rcut = args.rcut*angstrom
    alpha = args.alpha_scale/rcut
    gcut = args.gcut_scale*alpha
    return rcut, alpha, gcut


def parse_pbc(spbc):
    if len(spbc) != 3:
        raise ValueError('The pbc argument must consist of three characters, 0 or 1.')
    result = np.zeros(3, int)
    for i in xrange(3):
        if spbc[i] not in '01':
            raise ValueError('The pbc argument must consist of three characters, 0 or 1.')
        result[i] = int(spbc[i])
    return result


def parse_ugrid(sgrid, cell):
    '''Create a uniform integration grid based on the sgrid argument

       **Arguments:**

       sgrid
            A string obtained from the command line

       rvecs
            The 3D cell vectors
    '''
    try:
        spacing = float(sgrid)*angstrom
    except ValueError:
        spacing = None

    if spacing is None:
        fn_h5, path = sgrid.split(':')
        with h5.File(fn_h5, 'r') as f:
            return UniformGrid.from_hdf5(f[path], None)
    else:
        grid_rvecs = cell.rvecs.copy()
        shape = np.zeros(3, int)
        lengths, angles = cell.parameters
        for i in xrange(3):
            shape[i] = int(np.round(lengths[i]/spacing))
            grid_rvecs[i] /= shape[i]
        grid_cell = Cell(grid_rvecs)
        origin = np.zeros(3, float)
        pbc = np.ones(3, int)
        return UniformGrid(origin, grid_rvecs, shape, pbc)


def store_args(args, grp):
    '''Convert the command line arguments to hdf5 attributes'''
    grp.attrs['cmdline'] =  ' '.join(sys.argv)
    grp.attrs['pwd'] = os.getcwd()
    grp.attrs['datetime'] = datetime.datetime.now().isoformat()
    for key, val in vars(args).iteritems():
        if val is not None:
            grp.attrs['arg_%s' % key] = val


@contextlib.contextmanager
def safe_open_h5(*args, **kwargs):
    '''Try to open a file and 10 times wait a random time up to seconds between each attempt.

       **Arguments:** the same as those of the h5py.File constructor.

       **Optional keyword arguments:**

       wait
            The maximum number of seconds to wait between two attempts to open
            the file. [default=10]

       count
            The number of attempts to open the file.
    '''
    wait = kwargs.pop('wait', 10.0)
    counter = kwargs.pop('count', 10)
    while True:
        try:
            f = h5.File(*args, **kwargs)
            break
        except:
            counter -= 1
            if counter <= 0:
                raise
            time.sleep(np.random.uniform(0, wait))
            continue
    try:
        yield f
    finally:
        f.close()


def write_part_output(fn_h5, label, part, grp_name, names, args):
    # Store the results in an HDF5 file
    with safe_open_h5(fn_h5) as f:
        # Store system
        sys_grp = f.require_group('system')
        if 'cube_data' in part.system.props:
            del part.system.props['cube_data'] # drop potentially large array
        part.system.to_file(sys_grp)

        # Store results
        grp_part = f.require_group(label)
        if grp_name in grp_part:
            del grp_part[grp_name]
        grp = grp_part.create_group(grp_name)
        for name in names:
            dump_hdf5_low(grp, name, part[name])

        # Store command line arguments
        store_args(args, grp)

        if args.debug:
            # Store additional data for debugging
            if 'debug' in grp:
                del grp['debug']
            grp_debug = grp.create_group('debug')
            for debug_key in part.cache._store:
                debug_name = '_'.join(str(x) for x in debug_key)
                if debug_name not in names:
                    grp_debug[debug_name] = part.cache.load(*debug_key)

        if log.do_medium:
            log('Results written to %s:%s/%s' % (fn_h5, label, grp_name))
