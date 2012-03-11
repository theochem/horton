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
'''Input and output routines.

   All input routines begin with ``load_``. All output routines begin with
   ``dump_``.

   This module also contains a smart routine, ``load_system_args``, which makes
   it easier to load molecular data from various file formats. It uses to
   extension of the filename to figure out what the file format is.
'''


import numpy as np
from horton.units import angstrom
from horton.periodic import periodic


__all__ = [
    'load_system_args', 'load_geom_xyz', 'load_operators_g09', 'FCHKFile'
]


def load_system_args(filename):
    '''Load a molecular data from a file.

       **Argument:**

       filename
            The file to load the geometry from

       This routine uses the extension of the filename to determine the file
       format. It returns a dictionary with constructor arguments for the System
       class.
    '''
    if filename.endswith('.xyz'):
        coordinates, numbers = load_geom_xyz(filename)
        return {'coordinates': coordinates, 'numbers': numbers}
    else:
        raise ValueError('Unknown file format: %s' % filename)


def load_geom_xyz(filename):
    '''Load a molecular geometry from a .xyz file.

       **Argument:**

       filename
            The file to load the geometry from

       **Returns:** two arrays, coordinates and numbers that can be used as the
       two first arguments of the System constructor.
    '''
    f = file(filename)
    size = int(f.next())
    f.next()
    coordinates = np.empty((size, 3), float)
    numbers = np.empty(size, int)
    for i in xrange(size):
        words = f.next().split()
        numbers[i] = periodic[words[0]].number
        coordinates[i,0] = float(words[1])*angstrom
        coordinates[i,1] = float(words[2])*angstrom
        coordinates[i,2] = float(words[3])*angstrom
    f.close()
    return coordinates, numbers


def load_operators_g09(fn, lf):
    """Loads several one- and two-body operators from a Gaussian log file.

       **Arugment:**

       fn
            The filename of the Gaussian log file.

       lf
            A LinalgFactory instance.

       The following one-body operators are loaded if present: overlap, kinetic,
       nuclear attraction. The following two-body operator is loaded if present:
       electrostatic repulsion. In order to make all these matrices are present
       in the Gaussian log file, the following commands must be used in the
       Gaussian input file:

            scf(conventional) iop(3/33=5) extralinks=l316 iop(3/27=999)
    """


    with open(fn) as f:
        # First get the line with the number of basis functions
        for line in f:
            if line.startswith('    NBasis ='):
                nbasis = int(line[12:18])
                break

        # Then load the one- and two-body operators. This part is written such
        # that it does not make any assumptions about the order in which these
        # operators are printed.
        overlap = None
        kinetic = None
        nuclear_attraction = None
        electron_repulsion = None

        for line in f:
            if line == ' *** Overlap ***\n':
                overlap = _load_onebody_g09(f, nbasis, lf)
            elif line == ' *** Kinetic Energy ***\n':
                kinetic = _load_onebody_g09(f, nbasis, lf)
            elif line == ' ***** Potential Energy *****\n':
                nuclear_attraction = _load_onebody_g09(f, nbasis, lf)
            elif line == ' *** Dumping Two-Electron integrals ***\n':
                electron_repulsion = _load_twobody_g09(f, nbasis, lf)

        return overlap, kinetic, nuclear_attraction, electron_repulsion


def _load_onebody_g09(f, nbasis, lf):
    """Load a one-body operator from a Gaussian log file

       **Arguments:**

       f
            A file object for the Gaussian log file in read mode.

       nbasis
            The number of basis functions.

       lf
            A LinalgFactory instance.
    """
    result = lf.create_one_body(nbasis)
    block_counter = 0
    while block_counter < nbasis:
        # skip the header line
        f.next()
        # determine the number of rows in this part
        nrow = nbasis - block_counter
        for i in xrange(nrow):
            words = f.next().split()[1:]
            for j in xrange(len(words)):
                value = float(words[j].replace('D', 'E'))
                result.set_element(i+block_counter, j+block_counter, value)
        block_counter += 5
    return result


def _load_twobody_g09(f, nbasis, lf):
    """Load a two-body operator from a Gaussian log file

       **Arguments:**

       f
            A file object for the Gaussian log file in read mode.

       nbasis
            The number of basis functions.

       lf
            A LinalgFactory instance.
    """
    result = lf.create_two_body(nbasis)
    # Skip first six lines
    for i in xrange(6):
        f.next()
    # Start reading elements until a line is encountered that does not start
    # with ' I='
    while True:
        line = f.next()
        if not line.startswith(' I='):
            break
        #print line[3:7], line[9:13], line[15:19], line[21:25], line[28:].replace('D', 'E')
        i = int(line[3:7])-1
        j = int(line[9:13])-1
        k = int(line[15:19])-1
        l = int(line[21:25])-1
        value = float(line[29:].replace('D', 'E'))
        # Gaussian uses the chemists notation for the 4-center indexes. Horton
        # uses the physicists notation.
        result.set_element(i, k, j, l, value)
    return result


class FCHKFile(object):
    """Reader for Formatted checkpoint files

       After initialization, the data from the file is available in the fields
       dictionary. Also the following attributes are read from the file: title,
       command, lot (level of theory) and basis.
    """

    def __init__(self, filename, ignore_errors=False, field_labels=None):
        """
           **Arguments:**

           filename
                The formatted checkpoint file.

           **Optional arguments:**

           ignore_errors
                Try to read incorrectly formatted files without raising
                exceptions [default=False]

           field_labels
                When provided, only these fields are read from the formatted
                checkpoint file. (This can save a lot of time.)
        """
        self.filename = filename
        self.ignore_errors = ignore_errors
        try:
            self._read(filename, field_labels)
        except IOError:
            if ignore_errors:
                pass
            else:
                raise

    def _read(self, filename, field_labels=None):
        """Read all the requested fields"""
        # if fields is None, all fields are read
        def read_field(f):
            """Read a single field"""
            datatype = None
            while datatype is None:
                # find a sane header line
                line = f.readline()
                if line == "":
                    return False

                label = line[:43].strip()
                if field_labels is not None:
                    if len(field_labels) == 0:
                        return False
                    elif label not in field_labels:
                        return True
                    else:
                        field_labels.discard(label)
                line = line[43:]
                words = line.split()
                if len(words) == 0:
                    return True

                if words[0] == 'I':
                    datatype = np.int64
                    unreadable = 0
                elif words[0] == 'R':
                    datatype = float
                    unreadable = np.nan

            if len(words) == 2:
                try:
                    value = datatype(words[1])
                except ValueError:
                    return True
            elif len(words) == 3:
                if words[1] != "N=":
                    raise IOError("Unexpected line in formatted checkpoint file %s\n%s" % (filename, line[:-1]))
                length = int(words[2])
                value = np.zeros(length, datatype)
                counter = 0
                try:
                    while counter < length:
                        line = f.readline()
                        if line == "":
                            raise IOError("Unexpected end of formatted checkpoint file %s" % filename)
                        for word in line.split():
                            try:
                                value[counter] = datatype(word)
                            except (ValueError, OverflowError), e:
                                print 'WARNING: could not interpret word while reading %s: %s' % (word, self.filename)
                                if self.ignore_errors:
                                    value[counter] = unreadable
                                else:
                                    raise
                            counter += 1
                except ValueError:
                    return True
            else:
                raise IOError("Unexpected line in formatted checkpoint file %s\n%s" % (filename, line[:-1]))

            self.fields[label] = value
            return True

        self.fields = {}
        f = file(filename, 'r')
        self.title = f.readline()[:-1].strip()
        words = f.readline().split()
        if len(words) == 3:
            self.command, self.lot, self.basis = words
        elif len(words) == 2:
            self.command, self.lot = words
        else:
            raise IOError('The second line of the FCHK file should contain two or three words.')

        while read_field(f):
            pass

        f.close()
