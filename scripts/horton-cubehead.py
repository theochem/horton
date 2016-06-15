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


import argparse, numpy as np
from horton import IOData, __version__, angstrom


def parse_args():
    parser = argparse.ArgumentParser(prog='horton-cubehead.py',
        description='Write an ecomonic header for the cubegen program.')
    parser.add_argument('-V', '--version', action='version',
        version="%%(prog)s (HORTON version %s)" % __version__)

    parser.add_argument('structure',
        help='A file containing a molecular structure (of an isolated system)')
    parser.add_argument('output',
        help='The output filename to write the header to.')

    parser.add_argument('--spacing', default=0.2, type=float,
        help='The grid spacing in angstrom [default=%(default)s]')
    parser.add_argument('--margin', default=5.0, type=float,
        help='The margin between the box edges and the nuclei in angstrom '
             '[default=%(default)s]')

    return parser.parse_args()


def main():
    args = parse_args()
    margin = args.margin*angstrom
    spacing = args.spacing*angstrom

    mol = IOData.from_file(args.structure)
    # compute the shape tensor
    shape = np.dot(mol.coordinates.transpose(), mol.coordinates)
    # diagonalize to obtain the x, y and z directions.
    evals, evecs = np.linalg.eigh(shape)
    axes = evecs.transpose()*spacing

    # compute the origin and the number of repetitions along each axis.
    nrep = np.zeros(3, int)
    origin = np.zeros(3, float)
    for i in xrange(3):
        projc = np.dot(mol.coordinates, evecs[:,i])
        nrep[i] = np.ceil((projc.max() - projc.min() + 2*margin)/spacing)+1
        origin += 0.5*(projc.max() + projc.min() - (nrep[i]-1)*spacing)*evecs[:,i]

    with open(args.output, 'w') as f:
        # the header is written in Bohr, hence the -nrep[0]
        print >> f, '% 5i % 15.10f % 15.10f % 15.10f' % (0, origin[0], origin[1], origin[2])
        print >> f, '% 5i % 15.10f % 15.10f % 15.10f' % (-nrep[0], axes[0,0], axes[0,1], axes[0,2])
        print >> f, '% 5i % 15.10f % 15.10f % 15.10f' % (nrep[1], axes[1,0], axes[1,1], axes[1,2])
        print >> f, '% 5i % 15.10f % 15.10f % 15.10f' % (nrep[2], axes[2,0], axes[2,1], axes[2,2])


if __name__ == '__main__':
    main()
