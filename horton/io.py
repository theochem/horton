# -*- coding: utf-8 -*-
# Horton is a Density Functional Theory program.
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
'''Input and output routines.

   All input routines begin with ``load_``. All output routines begin with
   ``dump_``.

   This module also contains a smart routine, ``load_system_args``, which makes
   it easier to load molecular data from various file formats. It uses to
   extension of the filename to figure out what the file format is.
'''


import numpy as np, h5py as h5
from horton.units import angstrom
from horton.periodic import periodic


__all__ = [
    'load_system_args', 'load_geom_xyz', 'load_operators_g09', 'FCHKFile',
    'load_fchk',
]


def load_system_args(filename, lf):
    '''Load a molecular data from a file.

       **Argument:**

       filename
            The file to load the geometry from

       lf
            A LinalgFactory instance.

       This routine uses the extension of the filename to determine the file
       format. It returns a dictionary with constructor arguments for the System
       class.
    '''
    if isinstance(filename, h5.File) or filename.endswith('.h5'):
        return load_checkpoint(filename, lf)
    if filename.endswith('.xyz'):
        coordinates, numbers = load_geom_xyz(filename)
        return {'coordinates': coordinates, 'numbers': numbers}
    elif filename.endswith('.fchk'):
        coordinates, numbers, obasis, wfn, permutation, props, operators = load_fchk(filename, lf)
        return {'coordinates': coordinates, 'numbers': numbers, 'obasis': obasis,
                'wfn': wfn, 'permutation': permutation, 'props': props, 'operators': operators}
    elif filename.endswith('.log'):
        overlap, kinetic, nuclear_attraction, electronic_repulsion = load_operators_g09(filename, lf)
        operators = {}
        if overlap is not None:
            operators['olp'] = overlap
        if kinetic is not None:
            operators['kin'] = kinetic
        if nuclear_attraction is not None:
            operators['na'] = nuclear_attraction
        if electronic_repulsion is not None:
            operators['er'] = electronic_repulsion
        return {'operators': operators}
    elif filename.endswith('.mkl'):
        coordinates, numbers, obasis, wfn, signs = load_mkl(filename, lf)
        return {
            'coordinates': coordinates, 'numbers': numbers, 'obasis': obasis,
            'wfn': wfn, 'signs': signs,
        }
    elif filename.endswith('.molden.input'):
        coordinates, numbers, obasis, wfn, signs = load_molden(filename, lf)
        return {
            'coordinates': coordinates, 'numbers': numbers, 'obasis': obasis,
            'wfn': wfn, 'signs': signs,
        }
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
        # First get the line with the number of orbital basis functions
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
        electronic_repulsion = None

        for line in f:
            if line.startswith(' *** Overlap ***'):
                overlap = _load_onebody_g09(f, nbasis, lf)
            elif line.startswith(' *** Kinetic Energy ***'):
                kinetic = _load_onebody_g09(f, nbasis, lf)
            elif line.startswith(' ***** Potential Energy *****'):
                nuclear_attraction = _load_onebody_g09(f, nbasis, lf)
            elif line.startswith(' *** Dumping Two-Electron integrals ***'):
                electronic_repulsion = _load_twobody_g09(f, nbasis, lf)

        return overlap, kinetic, nuclear_attraction, electronic_repulsion


def _load_onebody_g09(f, nbasis, lf):
    """Load a one-body operator from a Gaussian log file

       **Arguments:**

       f
            A file object for the Gaussian log file in read mode.

       nbasis
            The number of orbital basis functions.

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
            The number of orbital basis functions.

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

    def __init__(self, filename, field_labels=None):
        """
           **Arguments:**

           filename
                The formatted checkpoint file.

           **Optional arguments:**

           field_labels
                When provided, only these fields are read from the formatted
                checkpoint file. (This can save a lot of time.)
        """
        self.filename = filename
        self._read(filename, set(field_labels))

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
                    datatype = int
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
            self.command, self.lot, self.obasis = words
        elif len(words) == 2:
            self.command, self.lot = words
        else:
            raise IOError('The second line of the FCHK file should contain two or three words.')

        while read_field(f):
            pass

        f.close()


def load_fchk(filename, lf):
    '''Load from a formatted checkpoint file.

       **Arguments:**

       filename
            The filename of the Gaussian formatted checkpoint file.

       lf
            A LinalgFactory instance.
    '''
    from horton.gbasis import GOBasis

    fchk = FCHKFile(filename, [
        "Number of electrons", "Number of independant functions",
        "Number of alpha electrons", "Number of beta electrons",
        "Atomic numbers", "Current cartesian coordinates",
        "Shell types", "Shell to atom map", "Shell to atom map",
        "Number of primitives per shell", "Primitive exponents",
        "Contraction coefficients", "P(S=P) Contraction coefficients",
        "Alpha Orbital Energies", "Alpha MO coefficients",
        "Beta Orbital Energies", "Beta MO coefficients",
        "Total Energy",
        'Total SCF Density', 'Spin SCF Density',
        'Total CC Density', 'Spin CC Density',
    ])

    # A) Load the geometry
    numbers = fchk.fields["Atomic numbers"]
    coordinates = fchk.fields["Current cartesian coordinates"].reshape(-1,3)

    # B) Load the orbital basis set
    shell_types = fchk.fields["Shell types"]
    shell_map = fchk.fields["Shell to atom map"] - 1
    nprims = fchk.fields["Number of primitives per shell"]
    alphas = fchk.fields["Primitive exponents"]
    ccoeffs_level1 = fchk.fields["Contraction coefficients"]
    ccoeffs_level2 = fchk.fields.get("P(S=P) Contraction coefficients")

    my_shell_types = []
    my_shell_map = []
    my_nprims = []
    my_alphas = []
    con_coeffs = []
    counter = 0
    for i, n in enumerate(nprims):
        if shell_types[i] == -1:
            # Special treatment for SP shell type
            my_shell_types.append(0)
            my_shell_types.append(1)
            my_shell_map.append(shell_map[i])
            my_shell_map.append(shell_map[i])
            my_nprims.append(nprims[i])
            my_nprims.append(nprims[i])
            my_alphas.append(alphas[counter:counter+n])
            my_alphas.append(alphas[counter:counter+n])
            con_coeffs.append(ccoeffs_level1[counter:counter+n])
            con_coeffs.append(ccoeffs_level2[counter:counter+n])
        else:
            my_shell_types.append(shell_types[i])
            my_shell_map.append(shell_map[i])
            my_nprims.append(nprims[i])
            my_alphas.append(alphas[counter:counter+n])
            con_coeffs.append(ccoeffs_level1[counter:counter+n])
        counter += n
    my_shell_types = np.array(my_shell_types)
    my_shell_map = np.array(my_shell_map)
    my_nprims = np.array(my_nprims)
    my_alphas = np.concatenate(my_alphas)
    con_coeffs = np.concatenate(con_coeffs)
    del shell_map
    del shell_types
    del nprims
    del alphas

    obasis = GOBasis(coordinates, my_shell_map, my_nprims, my_shell_types, my_alphas, con_coeffs)

    # permutation of the orbital basis functions
    permutation_rules = {
      -9: np.arange(19),
      -8: np.arange(17),
      -7: np.arange(15),
      -6: np.arange(13),
      -5: np.arange(11),
      -4: np.arange(9),
      -3: np.arange(7),
      -2: np.arange(5),
       0: np.array([0]),
       1: np.arange(3),
       2: np.array([0, 3, 4, 1, 5, 2]),
       3: np.array([0, 4, 5, 3, 9, 6, 1, 8, 7, 2]),
       4: np.arange(15)[::-1],
       5: np.arange(21)[::-1],
       6: np.arange(28)[::-1],
       7: np.arange(36)[::-1],
       8: np.arange(45)[::-1],
       9: np.arange(55)[::-1],
    }
    offset = 0
    permutation = []
    for shell_type in my_shell_types:
        permutation.extend(permutation_rules[shell_type]+len(permutation))
    permutation = np.array(permutation, dtype=int)

    # C) Load the wavefunction
    nbasis_indep = fchk.fields.get("Number of independant functions")
    if nbasis_indep is None:
        nbasis_indep = obasis.nbasis
    if 'Beta Orbital Energies' in fchk.fields:
        from horton.wfn import AufbauOccModel, OpenShellWFN
        nalpha = fchk.fields['Number of alpha electrons']
        nbeta = fchk.fields['Number of beta electrons']
        occ_model = AufbauOccModel(nalpha, nbeta)
        wfn = OpenShellWFN(occ_model, lf, obasis.nbasis, norb=nbasis_indep)
        exp_alpha = wfn.init_exp('alpha')
        exp_alpha.coeffs[:] = fchk.fields['Alpha MO coefficients'].reshape(nbasis_indep, obasis.nbasis).T
        exp_alpha.energies[:] = fchk.fields['Alpha Orbital Energies']
        exp_beta = wfn.init_exp('beta')
        exp_beta.coeffs[:] = fchk.fields['Beta MO coefficients'].reshape(nbasis_indep, obasis.nbasis).T
        exp_beta.energies[:] = fchk.fields['Beta Orbital Energies']
        occ_model.assign(exp_alpha, exp_beta)
    else:
        from horton.wfn import AufbauOccModel, ClosedShellWFN
        nelec = fchk.fields["Number of electrons"]
        assert nelec % 2 == 0
        occ_model = AufbauOccModel(nelec/2)
        wfn = ClosedShellWFN(occ_model, lf, obasis.nbasis, norb=nbasis_indep)
        exp_alpha = wfn.init_exp('alpha')
        exp_alpha.coeffs[:] = fchk.fields['Alpha MO coefficients'].reshape(nbasis_indep, obasis.nbasis).T
        exp_alpha.energies[:] = fchk.fields['Alpha Orbital Energies']
        occ_model.assign(exp_alpha)

    # D) Load properties
    props = {
        'energy': fchk.fields['Total Energy'],
    }

    # E) Load density matrices
    def load_dm(key, label):
        if label in fchk.fields:
            dm = lf.create_one_body(obasis.nbasis)
            start = 0
            for i in xrange(obasis.nbasis):
                stop = start+i+1
                dm._array[i,:i+1] = fchk.fields[label][start:stop]
                dm._array[:i+1,i] = fchk.fields[label][start:stop]
                start = stop
            operators[key] = dm

    # Note that the density matrices are not directly loaded into the
    # wavefunction objects. It is up to the user to decide which
    # density matrix will be used in the wavefunction, by making a few manual
    # assignments.
    operators = {}
    # TODO add more
    load_dm('scf_full', 'Total SCF Density')
    load_dm('scf_spin', 'Spin SCF Density')
    load_dm('mp2_full', 'Total MP2 Density')
    load_dm('mp2_spin', 'Spin MP2 Density')
    load_dm('cc_full', 'Total CC Density')
    load_dm('cc_spin', 'Spin CC Density')

    return coordinates, numbers, obasis, wfn, permutation, props, operators


def renorm_helper(con_coeff, alpha, shell_type):
    '''Fix an unnormalized contraction coefficient'''
    from horton.gbasis.cext import gob_cart_normalization
    if shell_type == 0:
        return con_coeff/gob_cart_normalization(alpha, np.array([0,0,0]))
    elif shell_type == 1:
        return con_coeff/gob_cart_normalization(alpha, np.array([1,0,0]))
    elif shell_type == -2:
        return con_coeff/gob_cart_normalization(alpha, np.array([1,1,0]))
    elif shell_type == -3:
        return con_coeff/gob_cart_normalization(alpha, np.array([1,1,1]))
    elif shell_type == -4:
        return con_coeff/gob_cart_normalization(alpha, np.array([2,1,1]))
    else:
        raise NotImplementedError('MKL Normalization conventions beyond G are not know. Please notify Toon.Verstraelen@UGent.be.')


def get_orca_signs(obasis):
    '''Correct for different sign conventions used in ORCA.'''
    sign_rules = {
      -4: np.array([1,1,1,1,1,-1,-1,-1,-1]),
      -3: np.array([1,1,1,1,1,-1,-1]),
      -2: np.array([1,1,1,1,1]),
       0: np.array([1]),
       1: np.array([1,1,1]),
    }
    signs = []
    for shell_type in obasis.shell_types:
        signs.extend(sign_rules[shell_type])
    return np.array(signs, dtype=int)


def load_molden(filename, lf):
    '''Load data from a molden input file (with ORCA sign conventions).

       **Arguments:**

       filename
            The filename of the molden input file.

       lf
            A LinalgFactory instance.
    '''

    def helper_coordinates(f):
        numbers = []
        coordinates = []
        while True:
            last_pos = f.tell()
            line = f.readline()
            if len(line) == 0:
                break
            words = line.split()
            if len(words) != 6:
                # Go back to previous line and stop
                f.seek(last_pos)
                break
            else:
                numbers.append(int(words[2]))
                coordinates.append([float(words[3]), float(words[4]), float(words[5])])
        numbers = np.array(numbers, int)
        coordinates = np.array(coordinates)
        return numbers, coordinates


    def helper_obasis(f, coordinates):
        from horton.gbasis.io import str_to_shell_types
        from horton.gbasis.gobasis import GOBasis
        shell_types = []
        shell_map = []
        nprims = []
        alphas = []
        con_coeffs = []

        icenter = 0
        in_atom = False
        in_shell = False
        while True:
            last_pos = f.tell()
            line = f.readline()
            if len(line) == 0:
                break
            words = line.split()
            if len(words) == 0:
                in_atom = False
                in_shell = False
            elif len(words) == 2 and not in_atom:
                icenter = int(words[0])-1
                in_atom = True
                in_shell = False
            elif len(words) == 3:
                in_shell = True
                shell_map.append(icenter)
                shell_type = str_to_shell_types(words[0])[0]
                if shell_type >= 2:
                    # always assume pure basis functions
                    shell_type *= -1
                shell_types.append(shell_type)
                nprims.append(int(words[1]))
            elif len(words) == 2 and in_atom:
                assert in_shell
                alpha = float(words[0])
                alphas.append(alpha)
                con_coeff = renorm_helper(float(words[1]), alpha, shell_type)
                con_coeffs.append(con_coeff)
            else:
                # done, go back one line
                f.seek(last_pos)
                break

        shell_map = np.array(shell_map)
        nprims = np.array(nprims)
        shell_types = np.array(shell_types)
        alphas = np.array(alphas)
        con_coeffs = np.array(con_coeffs)
        return GOBasis(coordinates, shell_map, nprims, shell_types, alphas, con_coeffs)


    def helper_coeffs(f, nbasis):
        coeff_alpha = []
        ener_alpha = []
        occ_alpha = []
        coeff_beta = []
        ener_beta = []
        occ_beta = []

        icoeff = -1
        while True:
            line = f.readline()
            if icoeff == nbasis:
                icoeff = -1
            if len(line) == 0:
                break
            elif icoeff == -1:
                # read 1a line
                if line != ' Sym=     1a\n':
                    raise IOError('Symmetry in wavefunctions is not supported.')
                # prepare array with orbital coefficients
                col = np.zeros((nbasis,1), float)
                icoeff = -2
            elif icoeff == -2:
                # read energy
                assert line.startswith(' Ene=')
                energy = float(line[5:])
                icoeff = -3
            elif icoeff == -3:
                # read expansion coefficients
                assert line.startswith(' Spin=')
                spin = line[6:].strip()
                if spin == 'Alpha':
                    coeff_alpha.append(col)
                    ener_alpha.append(energy)
                else:
                    coeff_beta.append(col)
                    ener_beta.append(energy)
                icoeff = -4
            elif icoeff == -4:
                assert line.startswith(' Occup=')
                occ = float(line[7:])
                if spin == 'Alpha':
                    occ_alpha.append(occ)
                else:
                    occ_beta.append(occ)
                icoeff = 0
            elif icoeff >= 0:
                words = line.split()
                col[icoeff] = float(words[1])
                icoeff+=1
                assert int(words[0]) == icoeff

        coeff_alpha = np.hstack(coeff_alpha)
        ener_alpha = np.array(ener_alpha)
        occ_alpha = np.array(occ_alpha)
        if len(coeff_beta) == 0:
            coeff_beta = None
            ener_beta = None
            occ_beta = None
        else:
            coeff_beta = np.hstack(coeff_beta)
            ener_beta = np.array(ener_beta)
            occ_beta = np.array(occ_beta)
        return (coeff_alpha, ener_alpha, occ_alpha), (coeff_beta, ener_beta, occ_beta)

    numbers = None
    coordinates = None
    obasis = None
    coeff_alpha = None
    ener_alpha = None
    occ_alpha = None
    coeff_beta = None
    ener_beta = None
    occ_beta = None
    with open(filename) as f:
        line = f.readline()
        if line != '[Molden Format]\n':
            raise IOError('Molden header not found')
        while True:
            line = f.readline()
            if len(line) == 0:
                break
            line = line.strip()
            if line == '[Atoms] AU':
                numbers, coordinates = helper_coordinates(f)
            elif line == '[GTO]':
                obasis = helper_obasis(f, coordinates)
            elif line == '[MO]':
                data_alpha, data_beta = helper_coeffs(f, obasis.nbasis)
                coeff_alpha, ener_alpha, occ_alpha = data_alpha
                coeff_beta, ener_beta, occ_beta = data_beta

    if coordinates is None:
        raise IOError('Coordinates not found in mkl file.')
    if obasis is None:
        raise IOError('Orbital basis not found in mkl file.')
    if coeff_alpha is None:
        raise IOError('Alpha orbitals not found in mkl file.')

    if coeff_beta is None:
        from horton.wfn import AufbauOccModel, ClosedShellWFN
        nalpha = int(np.round(occ_alpha.sum()))/2
        occ_model = AufbauOccModel(nalpha)
        wfn = ClosedShellWFN(occ_model, lf, obasis.nbasis, norb=coeff_alpha.shape[1])
        exp_alpha = wfn.init_exp('alpha')
        exp_alpha.coeffs[:] = coeff_alpha
        exp_alpha.energies[:] = ener_alpha
        exp_alpha.occupations[:] = occ_alpha/2
    else:
        from horton.wfn import AufbauOccModel, OpenShellWFN
        nalpha = int(np.round(occ_alpha.sum()))
        nbeta = int(np.round(occ_beta.sum()))
        assert coeff_alpha.shape == coeff_beta.shape
        assert ener_alpha.shape == ener_beta.shape
        assert occ_alpha.shape == occ_beta.shape
        occ_model = AufbauOccModel(nalpha, nbeta)
        wfn = OpenShellWFN(occ_model, lf, obasis.nbasis, norb=coeff_alpha.shape[1])
        exp_alpha = wfn.init_exp('alpha')
        exp_alpha.coeffs[:] = coeff_alpha
        exp_alpha.energies[:] = ener_alpha
        exp_alpha.occupations[:] = occ_alpha
        exp_beta = wfn.init_exp('beta')
        exp_beta.coeffs[:] = coeff_beta
        exp_beta.energies[:] = ener_beta
        exp_beta.occupations[:] = occ_beta

    signs = get_orca_signs(obasis)

    return coordinates, numbers, obasis, wfn, signs


def load_mkl(filename, lf):
    '''Load data from a Molekel file (with ORCA sign conventions).

       **Arguments:**

       filename
            The filename of the mkl file.

       lf
            A LinalgFactory instance.
    '''

    def helper_char_mult(f):
        return [int(word) for word in f.readline().split()]


    def helper_coordinates(f):
        numbers = []
        coordinates = []
        while True:
            line = f.readline()
            if len(line) == 0 or line.strip() == '$END':
                break
            words = line.split()
            numbers.append(int(words[0]))
            coordinates.append([float(words[1]), float(words[2]), float(words[3])])
        numbers = np.array(numbers, int)
        coordinates = np.array(coordinates)*angstrom
        return numbers, coordinates


    def helper_obasis(f, coordinates):
        from horton.gbasis.io import str_to_shell_types
        from horton.gbasis.gobasis import GOBasis
        shell_types = []
        shell_map = []
        nprims = []
        alphas = []
        con_coeffs = []

        center_counter = 0
        in_shell = False
        nprim = None
        while True:
            line = f.readline()
            lstrip = line.strip()
            if len(line) == 0 or lstrip == '$END':
                break
            if len(lstrip) == 0:
                continue
            if lstrip == '$$':
                center_counter += 1
                in_shell = False
            else:
                words = line.split()
                if len(words) == 2:
                    assert in_shell
                    alpha = float(words[0])
                    alphas.append(alpha)
                    con_coeff = renorm_helper(float(words[1]), alpha, shell_type)
                    con_coeffs.append(con_coeff)
                    nprim += 1
                else:
                    if nprim is not None:
                        nprims.append(nprim)
                    shell_map.append(center_counter)
                    shell_type = str_to_shell_types(words[1])[0]
                    if shell_type >= 2:
                        # always assume pure basis functions
                        shell_type *= -1
                    shell_types.append(shell_type)
                    in_shell = True
                    nprim = 0
        if nprim is not None:
            nprims.append(nprim)

        shell_map = np.array(shell_map)
        nprims = np.array(nprims)
        shell_types = np.array(shell_types)
        alphas = np.array(alphas)
        con_coeffs = np.array(con_coeffs)
        return GOBasis(coordinates, shell_map, nprims, shell_types, alphas, con_coeffs)


    def helper_coeffs(f, nbasis):
        coeffs = []
        energies = []

        in_orb = 0
        while True:
            line = f.readline()
            lstrip = line.strip()
            if len(line) == 0 or lstrip == '$END':
                break
            if in_orb == 0:
                # read a1g line
                words = lstrip.split()
                ncol = len(words)
                assert ncol > 0
                for word in words:
                    assert word == 'a1g'
                cols = [np.zeros((nbasis,1), float) for icol in xrange(ncol)]
                in_orb = 1
            elif in_orb == 1:
                # read energies
                words = lstrip.split()
                assert len(words) == ncol
                for word in words:
                    energies.append(float(word))
                in_orb = 2
                ibasis = 0
            elif in_orb == 2:
                # read expansion coefficients
                words = lstrip.split()
                assert len(words) == ncol
                for icol in xrange(ncol):
                    cols[icol][ibasis] = float(words[icol])
                ibasis += 1
                if ibasis == nbasis:
                    in_orb = 0
                    coeffs.extend(cols)

        return np.hstack(coeffs), np.array(energies)


    def helper_occ(f):
        occs = []
        while True:
            line = f.readline()
            lstrip = line.strip()
            if len(line) == 0 or lstrip == '$END':
                break
            for word in lstrip.split():
                occs.append(float(word))
        return np.array(occs)


    charge = None
    spinmult = None
    numbers = None
    coordinates = None
    obasis = None
    coeff_alpha = None
    ener_alpha = None
    occ_alpha = None
    coeff_beta = None
    ener_beta = None
    occ_beta = None
    with open(filename) as f:
        while True:
            line = f.readline()
            if len(line) == 0:
                break
            line = line.strip()
            if line == '$CHAR_MULT':
                charge, spinmult = helper_char_mult(f)
            elif line == '$COORD':
                numbers, coordinates = helper_coordinates(f)
            elif line == '$BASIS':
                obasis = helper_obasis(f, coordinates)
            elif line == '$COEFF_ALPHA':
                coeff_alpha, ener_alpha = helper_coeffs(f, obasis.nbasis)
            elif line == '$OCC_ALPHA':
                occ_alpha = helper_occ(f)
            elif line == '$COEFF_BETA':
                coeff_beta, ener_beta = helper_coeffs(f, obasis.nbasis)
            elif line == '$OCC_BETA':
                occ_beta = helper_occ(f)

    if charge is None:
        raise IOError('Charge and multiplicity not found in mkl file.')
    if coordinates is None:
        raise IOError('Coordinates not found in mkl file.')
    if obasis is None:
        raise IOError('Orbital basis not found in mkl file.')
    if coeff_alpha is None:
        raise IOError('Alpha orbitals not found in mkl file.')
    if occ_alpha is None:
        raise IOError('Alpha occupation numbers not found in mkl file.')

    nelec = numbers.sum() - charge
    if coeff_beta is None:
        from horton.wfn import AufbauOccModel, ClosedShellWFN
        assert nelec % 2 == 0
        assert abs(occ_alpha.sum() - nelec) < 1e-7
        occ_model = AufbauOccModel(nelec/2)
        wfn = ClosedShellWFN(occ_model, lf, obasis.nbasis, norb=coeff_alpha.shape[1])
        exp_alpha = wfn.init_exp('alpha')
        exp_alpha.coeffs[:] = coeff_alpha
        exp_alpha.energies[:] = ener_alpha
        exp_alpha.occupations[:] = occ_alpha/2
    else:
        if occ_beta is None:
            raise IOError('Beta occupation numbers not found in mkl file while beta orbitals were present.')
        from horton.wfn import AufbauOccModel, OpenShellWFN
        nalpha = int(np.round(occ_alpha.sum()))
        nbeta = int(np.round(occ_beta.sum()))
        assert nelec == nalpha+nbeta
        assert coeff_alpha.shape == coeff_beta.shape
        assert ener_alpha.shape == ener_beta.shape
        assert occ_alpha.shape == occ_beta.shape
        occ_model = AufbauOccModel(nalpha, nbeta)
        wfn = OpenShellWFN(occ_model, lf, obasis.nbasis, norb=coeff_alpha.shape[1])
        exp_alpha = wfn.init_exp('alpha')
        exp_alpha.coeffs[:] = coeff_alpha
        exp_alpha.energies[:] = ener_alpha
        exp_alpha.occupations[:] = occ_alpha
        exp_beta = wfn.init_exp('beta')
        exp_beta.coeffs[:] = coeff_beta
        exp_beta.energies[:] = ener_beta
        exp_beta.occupations[:] = occ_beta

        exp_beta.occupations[:] = occ_beta

    signs = get_orca_signs(obasis)

    return coordinates, numbers, obasis, wfn, signs


def load_checkpoint(filename, lf):
    """Load constructor arguments from a Horton checkpoint file

       **Arguments:**

       filename
            This is the file name of an HDF5 Horton checkpoint file. It may also
            be an open h5.File object.

       lf
            A LinalgFactory instance.
    """
    result = {}
    if isinstance(filename, basestring):
        chk = h5.File(filename,'r')
        do_close = True
    else:
        chk = filename
        do_close = False
    from horton.checkpoint import register
    for field in register.itervalues():
        att = field.read(chk, lf)
        d = result
        name = field.att_name
        if field.key is not None:
            d = result.setdefault(field.att_name, {})
            name = field.key
        if att is not None:
            d[name] = att
    if do_close:
        chk.close()
    result['chk'] = filename
    return result
