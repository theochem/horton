#!/usr/bin/env python

from glob import glob
import os
from horton import periodic

elements = set([])
grids = set([])
npoints = {}

for fn_txt in glob('../data/grids/*.txt'):
    grid = os.path.basename(fn_txt)[:-4]
    grids.add(grid)
    with open(fn_txt) as f:
        for line in f:
            if line.startswith('# Number of points:'):
                npoint = int(line.split()[-1])
            words = line.split()
            if line[0] != '#' and len(words) == 2:
                assert isinstance(npoint, int)
                npoints[int(words[0]), int(words[1]), grid] = npoint
                npoint = None
                elements.add((int(words[0]), int(words[1])))


grids = sorted(grids)
alts = ['coarse', 'medium', 'fine', 'veryfine', 'ultrafine', 'insane']
elements = sorted(elements)

with open('grids.rst.inc', 'w') as f:
    def print_markers():
        print >> f, '======= === ====', ' '.join(['=========']*len(grids))


    print_markers()
    print >> f, 'Element   Z Zeff', ' '.join([grid.ljust(9) for grid in grids])
    print >> f, '\         \ \   ', ' '.join([alt.ljust(9) for alt in alts])
    print_markers()

    for z, zeff in elements:
        print >> f, '%7s %3i %4i' % (periodic[z].symbol, z, zeff),
        for grid in grids:
            npoint = npoints.get((z, zeff, grid))
            if npoint is None:
                print >> f, '--       ',
            else:
                print >> f, '%9i' % npoint,
        print >> f

    print_markers()
