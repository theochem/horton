#!/usr/bin/env python
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


keep_fields = set([
    "Charge", "Multiplicity", "Number of electrons",
    "Number of alpha electrons", "Number of beta electrons",
    "Number of basis functions", "Number of independent functions",
    "Atomic numbers", "Nuclear charges", "Current cartesian coordinates",
    "Shell types", "Number of primitives per shell", "Shell to atom map",
    "Primitive exponents", "Contraction coefficients",
    "P(S=P) Contraction coefficients", "Coordinates of each shell",
    "Virial Ratio", "Total Energy", "Alpha Orbital Energies",
    "Beta Orbital Energies", "Alpha MO coefficients", "Beta MO coefficients",
    'Total SCF Density', 'Spin SCF Density',
    'Total MP2 Density', 'Spin MP2 Density',
    'Total MP3 Density', 'Spin MP3 Density',
    'Total CC Density', 'Spin CC Density',
    'Total CI Density', 'Spin CI Density',
    "Mulliken Charges", "Dipole Moment", "Quadrupole Moment",
    "Real atomic weights",
    "Cartesian Gradient", "Cartesian Force Constants",
])

def strip(fn):
    # a list with all lines that we'd like to keep
    lines = []
    # load the good lines
    f = file(fn, "r")
    busy = False
    keep = False
    for line in f:
        line = line.rstrip()
        if len(lines) < 2:
            # keep the first two lines
            lines.append(line)
        else:
            if busy:
                if not line.startswith(" "):
                    keep = False
                    busy = False
                if keep:
                    lines.append(line)
            if not busy:
                busy = line[47:49] == "N="
                field = line[:43].strip()
                if field in keep_fields:
                    lines.append(line)
                    keep = True
                else:
                    keep = False
    f.close()
    # print stuff back into the same file
    f = file(fn, "w")
    for line in lines:
        print >> f, line
    f.close()


if __name__ == "__main__":
    import sys
    fns_fchk = sys.argv[1:]
    for fn in fns_fchk:
        if fn.endswith(".fchk"):
            print "Stripping", fn
            strip(fn)
        else:
            print "Skipping", fn, "(wrong extension)"
