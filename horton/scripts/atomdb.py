# -*- coding: utf-8 -*-
# HORTON: Helpful Open-source Research TOol for N-fermion systems.
# Copyright (C) 2011-2016 The HORTON Development Team
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
"""Code used by ``horton-atomdb.py``"""


from glob import glob
import os
import re
import stat
from string import Template as BaseTemplate

import numpy as np
import matplotlib.pyplot as pt

from horton.io.iodata import IOData
from horton.log import log
from horton.periodic import periodic
from horton.scripts.common import iter_elements
from horton.units import angstrom


__all__ = [
    'iter_mults', 'iter_states', 'plot_atoms',
    'Template', 'EnergyTable', 'atom_programs',
]


# Presets for spin multiplicites. The first element is according to Hund's rule.
# Following elements are reasonable.
mult_presets = {
    1: [2],
    2: [1, 3],
    3: [2, 4],
    4: [1, 3],
    5: [2, 4],
    6: [3, 5, 1],
    7: [4, 2],
    8: [3, 1],
    9: [2],
    10: [1],
    11: [2],
    12: [1, 3],
    13: [2, 4],
    14: [3, 5, 1],
    15: [4, 2],
    16: [3, 1],
    17: [2],
    18: [1],
    19: [2],
    20: [1, 3],
    21: [2, 4],
    22: [3, 5, 1],
    23: [4, 6, 2],
    24: [7, 5, 3, 1],
    25: [6, 4, 2],
    26: [5, 3, 1],
    27: [4, 2],
    28: [3, 1],
    29: [2],
    30: [1],
    31: [2, 4, 6],
    32: [3, 1, 5],
    33: [4, 2],
    34: [3, 1],
    35: [2],
    36: [1],
    37: [2],
    38: [1, 3],
    39: [2, 4],
    40: [3, 1, 5],
    41: [6, 4, 2],
    42: [7, 5, 3, 1],
    43: [6, 4, 2],
    44: [5, 3, 1],
    45: [4, 2],
    46: [1, 3],
    47: [2],
    48: [1],
    49: [2, 4],
    50: [3, 1, 5],
    51: [4, 2],
    52: [3, 1],
    53: [2],
    54: [1],
    55: [2],
    56: [1, 3],
    57: [2, 4],
    58: [3, 1, 5],
    59: [4, 2, 6],
    60: [5, 1, 3, 7],
    61: [6, 2, 4, 8],
    62: [7, 1, 3, 5, 9],
    63: [8, 2, 4, 6, 10],
    64: [9, 1, 3, 5, 7, 11],
    65: [6, 2, 4, 8, 10, 12],
    66: [5, 1, 3, 7, 9, 11],
    67: [4, 2, 6, 8, 10],
    68: [3, 1, 5, 7, 9],
    69: [2, 4, 6, 8],
    70: [1, 2, 5, 7],
    71: [2, 4, 6],
    72: [3, 5, 1],
    73: [4, 2, 6],
    74: [5, 1, 3, 7],
    75: [6, 2, 4, 8],
    76: [5, 1, 3],
    77: [4, 2],
    78: [3, 1],
    79: [2],
    80: [1],
    81: [2, 4],
    82: [3, 1, 5],
    83: [4, 2],
    84: [3, 1],
    85: [2],
    86: [1],
}


def iter_mults(nel, hund):
    """Iterate over atomic spin multiplicites for the given number of electrons

       **Arguments:**

       nel
            The number of electrons (1-56)

       hund
            When set to True, only one spin multiplicity is returned. Otherwise
            several reasonable spin multiplicities are given.
    """

    if hund:
        yield mult_presets[nel][0]
    else:
        for mult in mult_presets[nel]:
            yield mult


def iter_states(elements, max_cation, max_anion, hund):
    """Iterate over all requested atomic states

       **Arguments:**

       elements
            A string that is suitable for ``iter_elements``

       template
            An instance of the ``atomdb.Template`` class

       max_cation
            The limit for the most positive cation

       max_anion
            The limit for the most negative anion

       hund
            Flag to adhere to hund's rule for the spin multiplicities.
    """
    for number in iter_elements(elements):
        # Loop over all charge states for this element
        for charge in xrange(-max_anion, max_cation+1):
            nel = number - charge
            if nel <= 0:
                continue
            # loop over multiplicities
            for mult in iter_mults(nel, hund):
                yield number, charge, mult


def plot_atoms(proatomdb, dn='.'):
    """Make PNG figures for all atoms in a pro-atom database.

    Warning: this script writes a bunch of PNG files!

    Parameters
    ----------
    proatomdb : horton.part.proatomdb.ProAtomDB
                A database of pro-atoms.
    dn : str
         Directory where the PNG files will be written. Local directory if not given.
    """
    def get_color(index):
        """Return a nice color for a given index."""
        colors = ["#FF0000", "#FFAA00", "#00AA00", "#00AAFF", "#0000FF", "#FF00FF", "#777777"]
        return colors[index % len(colors)]

    lss = {True: '-', False: ':'}
    for number in proatomdb.get_numbers():
        r = proatomdb.get_rgrid(number).radii
        symbol = periodic[number].symbol
        charges = proatomdb.get_charges(number)
        suffix = '%03i_%s' % (number, symbol.lower().rjust(2, '_'))

        # The density (rho)
        pt.clf()
        for i, charge in enumerate(charges):
            record = proatomdb.get_record(number, charge)
            y = record.rho
            ls = lss[record.safe]
            color = get_color(i)
            label = 'q=%+i' % charge
            pt.semilogy(r/angstrom, y, lw=2, ls=ls, label=label, color=color)
        pt.xlim(0, 3)
        pt.ylim(ymin=1e-5)
        pt.xlabel('Distance from the nucleus [A]')
        pt.ylabel('Spherically averaged density [Bohr**-3]')
        pt.title('Proatoms for element %s (%i)' % (symbol, number))
        pt.legend(loc=0)
        fn_png = '%s/dens_%s.png' % (dn, suffix)
        pt.savefig(fn_png)
        if log.do_medium:
            log('Written', fn_png)

        # 4*pi*r**2*rho
        pt.clf()
        for i, charge in enumerate(charges):
            record = proatomdb.get_record(number, charge)
            y = record.rho
            ls = lss[record.safe]
            color = get_color(i)
            label = 'q=%+i' % charge
            pt.plot(r/angstrom, 4*np.pi*r**2*y, lw=2, ls=ls, label=label, color=color)
        pt.xlim(0, 3)
        pt.ylim(ymin=0.0)
        pt.xlabel('Distance from the nucleus [A]')
        pt.ylabel('4*pi*r**2*density [Bohr**-1]')
        pt.title('Proatoms for element %s (%i)' % (symbol, number))
        pt.legend(loc=0)
        fn_png = '%s/rdens_%s.png' % (dn, suffix)
        pt.savefig(fn_png)
        if log.do_medium:
            log('Written', fn_png)

        fukui_data = []
        if number - charges[0] == 1:
            record0 = proatomdb.get_record(number, charges[0])
            fukui_data.append((record0.rho, record0.safe, '%+i' % charges[0]))
        for i, charge in enumerate(charges[1:]):
            record0 = proatomdb.get_record(number, charge)
            record1 = proatomdb.get_record(number, charges[i])
            fukui_data.append((
                record0.rho - record1.rho,
                record0.safe and record1.safe,
                '%+i-%+i' % (charge, charges[i])
            ))

        # The Fukui functions
        pt.clf()
        for i, (f, safe, label) in enumerate(fukui_data):
            ls = lss[safe]
            color = get_color(i)
            pt.semilogy(r/angstrom, f, lw=2, ls=ls, label=label, color=color, alpha=1.0)
            pt.semilogy(r/angstrom, -f, lw=2, ls=ls, color=color, alpha=0.2)
        pt.xlim(0, 3)
        pt.ylim(ymin=1e-5)
        pt.xlabel('Distance from the nucleus [A]')
        pt.ylabel('Fukui function [Bohr**-3]')
        pt.title('Proatoms for element %s (%i)' % (symbol, number))
        pt.legend(loc=0)
        fn_png = '%s/fukui_%s.png' % (dn, suffix)
        pt.savefig(fn_png)
        if log.do_medium:
            log('Written', fn_png)

        # 4*pi*r**2*Fukui
        pt.clf()
        for i, (f, safe, label) in enumerate(fukui_data):
            ls = lss[safe]
            color = get_color(i)
            pt.plot(r/angstrom, 4*np.pi*r**2*f, lw=2, ls=ls, label=label, color=color)
        pt.xlim(0, 3)
        pt.xlabel('Distance from the nucleus [A]')
        pt.ylabel('4*pi*r**2*Fukui [Bohr**-1]')
        pt.title('Proatoms for element %s (%i)' % (symbol, number))
        pt.legend(loc=0)
        fn_png = '%s/rfukui_%s.png' % (dn, suffix)
        pt.savefig(fn_png)
        if log.do_medium:
            log('Written', fn_png)


class Template(BaseTemplate):
    """A template with modifications to support inclusion of other files."""
    idpattern = r'[_a-z0-9.:-]+'

    def __init__(self, *args, **kwargs):
        BaseTemplate.__init__(self, *args, **kwargs)
        self._init_include_names()
        self._load_includes()

    def _init_include_names(self):
        """Return a list of include variables

           The include variables in the template are variables of the form
           ${file:name} or ${line:name}. This routine lists all the names
           encountered. Duplicates are eliminated.
        """
        pattern = '%s{(?P<braced>%s)}' % (re.escape(self.delimiter), self.idpattern)
        file_names = set([])
        line_names = set([])
        for mo in re.finditer(pattern, self.template):
            braced = mo.group('braced')
            if braced is not None and braced.startswith('file:'):
                file_names.add(braced[5:])
            if braced is not None and braced.startswith('line:'):
                line_names.add(braced[5:])
        self.file_names = list(file_names)
        self.line_names = list(line_names)

    def _load_includes(self):
        """Load included files for a given element number"""
        self.includes = []
        # Load files
        for name in self.file_names:
            records = []
            for fn in sorted(glob('%s.[0-9][0-9][0-9]_[0-9][0-9][0-9]_[0-9][0-9]' % name)):
                with open(fn) as f:
                    s = f.read()
                    # chop of one final newline if present (mostly the case)
                    if s[-1] == '\n':
                        s = s[:-1]
                number = int(fn[-10:-7])
                pop = int(fn[-6:-3])
                mult = int(fn[-2:])
                records.append((number, pop, mult, s))
            self.includes.append((name, 'file', records))
        # Load lines
        for name in self.line_names:
            with open(name) as f:
                records = []
                for line in f:
                    # ignore empty lines
                    if len(line.strip()) == 0:
                        continue
                    number = int(line[:3])
                    assert line[3] == '_'
                    pop = int(line[4:7])
                    assert line[7] == '_'
                    mult = int(line[8:10])
                    assert line[10] == ' '
                    s = line[11:-1]
                    records.append((number, pop, mult, s))
            self.includes.append((name, 'line', records))

    def _log_includes(self):
        # log the include names
        if len(self.file_names) + len(self.line_names) > 0 and log.do_medium:
            log('The following includes were detected in the template:')
            for name, kind, records in self.includes:
                log('   %s (%s)' % (name, kind))
                for n, p, m, s in self.includes[name]:
                    log('      %03i_%03i_%02i' % (n, p, m))

    def get_subs(self, number, pop, mult):
        subs = {}
        for name, kind, records in self.includes:
            found_s = None
            for n, p, m, s in records:
                if ((n==0) or (number==n)) and ((p==0) or (pop==p)) and ((m==0) or (mult==m)):
                    # match
                    found_s = s
                    break
            if found_s is None:
                raise KeyError('No matching include found for \'%s\' (%03i_%03i_%02i)' % (name, number, pop, mult))
            subs['%s:%s' % (kind, name)] = s
        return subs


class EnergyTable(object):
    def __init__(self):
        self.all = {}

    def add(self, number, pop, energy):
        cases = self.all.setdefault(number, {})
        cases[pop] = energy

    def log(self):
        log(' Nr Pop Chg              Energy          Ionization            Affinity')
        log.hline()
        for number, cases in sorted(self.all.iteritems()):
            for pop, energy in sorted(cases.iteritems()):
                energy_prev = cases.get(pop-1)
                if energy_prev is None:
                    ip_str = ''
                else:
                    ip_str = '% 18.10f' % (energy_prev - energy)

                energy_next = cases.get(pop+1)
                if energy_next is None:
                    ea_str = ''
                else:
                    ea_str = '% 18.10f' % (energy - energy_next)

                log('%3i %3i %+3i  % 18.10f  %18s  %18s' % (
                    number, pop, number-pop, energy, ip_str, ea_str
                ))
            log.blank()


class AtomProgram(object):
    name = None
    run_script = None

    def write_input(self, number, charge, mult, template, do_overwrite):
        # Directory stuff
        nel = number - charge
        dn_mult = '%03i_%s_%03i_q%+03i/mult%02i' % (
            number, periodic[number].symbol.lower().rjust(2, '_'), nel, charge, mult)

        # Figure out if we want to write
        fn_inp = '%s/atom.in' % dn_mult
        exists = os.path.isfile(fn_inp)
        do_write = not exists or do_overwrite

        if do_write:
            try:
                subs = template.get_subs(number, nel, mult)
            except KeyError:
                if log.do_warning:
                    log.warn('Could not find all subs for %03i.%03i.%03i. Skipping.' % (number, nel, mult))
                return dn_mult, False
            if not os.path.isdir(dn_mult):
                os.makedirs(dn_mult)
            with open(fn_inp, 'w') as f:
                f.write(template.substitute(
                    subs,
                    charge=str(charge),
                    mult=str(mult),
                    number=str(number),
                    element=periodic[number].symbol,
                ))
            if log.do_medium:
                if exists:
                    log('Overwritten:      ', fn_inp)
                else:
                    log('Written new:      ', fn_inp)
        elif log.do_medium:
            log('Not overwriting:  ', fn_inp)

        return dn_mult, do_write

    def write_run_script(self):
        # write the script
        fn_script = 'run_%s.sh' % self.name
        exists = os.path.isfile(fn_script)
        if not exists:
            with open(fn_script, 'w') as f:
                print >> f, self.run_script
            log('Written new:      ', fn_script)
        else:
            log('Not overwriting:  ', fn_script)
        # make the script executable
        os.chmod(fn_script, stat.S_IXUSR | os.stat(fn_script).st_mode)


    def _get_energy(self, mol, dn_mult):
        return mol.energy

    def load_atom(self, dn_mult, ext):
        fn = '%s/atom.%s' % (dn_mult, ext)
        if not os.path.isfile(fn):
            return None, None

        try:
            mol = IOData.from_file(fn)
        except:
            return None, None
        mol.energy = self._get_energy(mol, dn_mult)
        return mol, mol.energy


run_gaussian_script = """\
#!/bin/bash

# make sure %(name)s and formchk are available before running this script.

MISSING=0
if ! which %(name)s &>/dev/null; then echo "%(name)s binary not found."; MISSING=1; fi
if ! which formchk &>/dev/null; then echo "formchk binary not found."; MISSING=1; fi
if [ $MISSING -eq 1 ]; then echo "The required programs are not present on your system. Giving up."; exit -1; fi

function do_atom {
    echo "Computing in ${1}"
    cd ${1}
    if [ -e atom.out ]; then
        echo "Output file present in ${1}, not recomputing."
    else
        %(name)s atom.in > atom.out
        RETCODE=$?
        if [ $RETCODE == 0 ]; then
            formchk atom.chk atom.fchk
            rm -f atom.out.failed
        else
            # Rename the output of the failed job such that it gets recomputed
            # when the run script is executed again.
            mv atom.out atom.out.failed
        fi
        rm atom.chk
    fi
    cd -
}

for ATOMDIR in [01][0-9][0-9]_*_[01][0-9][0-9]_q[-+][0-9][0-9]/mult[0-9][0-9]; do
    do_atom ${ATOMDIR}
done
"""


class G09AtomProgram(AtomProgram):
    name = 'g09'
    run_script = run_gaussian_script % {'name': 'g09'}

    def write_input(self, number, charge, mult, template, do_overwrite):
        if '%chk=atom.chk\n' not in template.template:
            raise ValueError('The template must contain a line \'%chk=atom.chk\'')
        return AtomProgram.write_input(self, number, charge, mult, template, do_overwrite)

    def load_atom(self, dn_mult):
        return AtomProgram.load_atom(self, dn_mult, 'fchk')


class G03AtomProgram(G09AtomProgram):
    name = 'g03'
    run_script = run_gaussian_script % {'name': 'g03'}


run_orca_script = """\
#!/bin/bash

# make sure orca and orca2mkl are available before running this script.

MISSING=0
if ! which orca &>/dev/null; then echo "orca binary not found."; MISSING=1; fi
if ! which orca_2mkl &>/dev/null; then echo "orca_2mkl binary not found."; MISSING=1; fi
if [ $MISSING -eq 1 ]; then echo "The required programs are not present on your system. Giving up."; exit -1; fi

function do_atom {
    echo "Computing in ${1}"
    cd ${1}
    if [ -e atom.out ]; then
        echo "Output file present in ${1}, not recomputing."
    else
        orca atom.in > atom.out
        RETCODE=$?
        if [ $RETCODE == 0 ]; then
            orca_2mkl atom -molden
            rm -f atom.out.failed
        else
            # Rename the output of the failed job such that it gets recomputed
            # when the run script is executed again.
            mv atom.out atom.out.failed
        fi
    fi
    cd -
}

for ATOMDIR in [01][0-9][0-9]_*_[01][0-9][0-9]_q[-+][0-9][0-9]/mult[0-9][0-9]; do
    do_atom ${ATOMDIR}
done
"""


class OrcaAtomProgram(AtomProgram):
    name = 'orca'
    run_script = run_orca_script

    def _get_energy(self, mol, dn_mult):
        with open('%s/atom.out' % dn_mult) as f:
            for line in f:
                if line.startswith('Total Energy       :'):
                    return float(line[25:43])

    def load_atom(self, dn_mult):
        return AtomProgram.load_atom(self, dn_mult, 'molden.input')


run_cp2k_script = """\
#!/bin/bash

# Note: if you want to use an mpi-parallel CP2K binary, uncomment the following
# line and fill in the right binary and mpirun script:
#CP2K_BIN="mpirun -n4 cp2k.popt"

# Find a non-mpi CP2K binary if needed.
if [ -z "$CP2K_BIN" ]; then
    # Find all potential non-mpi CP2K binaries in the $PATH
    CP2K_BINS=$(find ${PATH//:/ } -name "cp2k.s*")

    # Check for any known non-mpi cp2k binary name in order of preference:
    for KNOWN_CP2K in cp2k.ssmp cp2k.sopt cp2k.sdbg; do
        for TMP in ${CP2K_BINS}; do
            if [ $(basename $TMP) == ${KNOWN_CP2K} ]; then
                CP2K_BIN=$TMP
                break
            fi
        done
        if [ -n $CP2K_BIN ]; then break; fi
    done

    MISSING=0
    if [ -z $CP2K_BIN ]; then echo "No non-mpi CP2K binary found."; MISSING=1; fi
    if [ $MISSING -eq 1 ]; then echo "The required programs are not present on your system. Giving up."; exit -1; fi
fi

echo "Using the following CP2K binary: $CP2K_BIN"

function do_atom {
    echo "Computing in ${1}"
    cd ${1}
    if [ -e atom.cp2k.out ]; then
        echo "Output file present in ${1}, not recomputing."
    else
        $CP2K_BIN atom.in > atom.cp2k.out
        RETCODE=$?
        if [ $RETCODE == 0 ]; then
            rm -f atom.cp2k.out.failed
        else
            # Rename the output of the failed job such that it gets recomputed
            # when the run script is executed again.
            mv atom.cp2k.out atom.cp2k.out.failed
        fi
    fi
    cd -
}

for ATOMDIR in [01][0-9][0-9]_*_[01][0-9][0-9]_q[-+][0-9][0-9]/mult[0-9][0-9]; do
    do_atom ${ATOMDIR}
done
"""

class CP2KAtomProgram(AtomProgram):
    name = 'cp2k'
    run_script = run_cp2k_script

    def write_input(self, number, charge, mult, template, do_overwrite):
        if '&ATOM' not in template.template:
            raise ValueError('The template must be a CP2K atom input. (\'&ATOM\' not found.)')
        return AtomProgram.write_input(self, number, charge, mult, template, do_overwrite)

    def load_atom(self, dn_mult):
        return AtomProgram.load_atom(self, dn_mult, 'cp2k.out')


run_psi4_script = """\
#!/bin/bash

# make sure psi4 is available before running this script.

MISSING=0
if ! which psi4 &>/dev/null; then echo "psi4 binary not found."; MISSING=1; fi
if [ $MISSING -eq 1 ]; then echo "The required programs are not present on your system. Giving up."; exit -1; fi

function do_atom {
    echo "Computing in ${1}"
    cd ${1}
    if [ -e atom.out ]; then
        echo "Output file present in ${1}, not recomputing."
    else
        psi4 atom.in
        RETCODE=$?
        if [ $RETCODE == 0 ]; then
            rm -f atom.out.failed
        else
            # Rename the output of the failed job such that it gets recomputed
            # when the run script is executed again.
            mv atom.out atom.out.failed
        fi
    fi
    cd -
}

for ATOMDIR in [01][0-9][0-9]_*_[01][0-9][0-9]_q[-+][0-9][0-9]/mult[0-9][0-9]; do
    do_atom ${ATOMDIR}
done
"""


class Psi4AtomProgram(AtomProgram):
    name = 'psi4'
    run_script = run_psi4_script

    def _get_energy(self, mol, dn_mult):
        with open('%s/atom.out' % dn_mult) as f:
            for line in f:
                if 'Final Energy' in line:
                    return float(line.split()[-1])

    def write_input(self, number, charge, mult, template, do_overwrite):
        found = False
        for line in template.template.split('\n'):
            words = line.lower().split()
            if 'molden_write' in words and 'true' in words:
                found = True
                break
        if not found:
            raise ValueError('The template must contain a line with \'molden_write true\'.')
        return AtomProgram.write_input(self, number, charge, mult, template, do_overwrite)

    def load_atom(self, dn_mult):
        return AtomProgram.load_atom(self, dn_mult, 'default.molden')


atom_programs = {}
for APC in globals().values():
    if isinstance(APC, type) and issubclass(APC, AtomProgram) and not APC is AtomProgram:
        atom_programs[APC.name] = APC()
