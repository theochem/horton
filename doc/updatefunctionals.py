#!/usr/bin/env python

libxc_version = '2.0.3'

# find all the functional keys by processing funcs_key.c

keys = []
with open('../depends/libxc-%s/src/funcs_key.c' % libxc_version) as f:
    for line in f:
        if line.startswith('{'):
            words = line.strip()[1:-3].split(',')
            key = words[0][1:-1]
            if len(key) > 0 and 'mgga' not in key:
                keys.append(key)

# make a rst table of all functionals

from horton import LibXCWrapper
from cStringIO import StringIO
from common import write_if_changed

s = StringIO()

print >> s, '.. _ref_functionals:'
print >> s
print >> s, 'LibXC Functionals'
print >> s, '#################'
print >> s
print >> s, 'The following functionals are available in Horton through `LibXC <http://www.tddft.org/programs/octopus/wiki/index.php/Libxc>`_ %s. [marques2012]_' % libxc_version
print >> s, '(The MGGA functionals are not supported yet in Horton.)'
print >> s
for key in keys:
    print >> s, '**%s**' % key
    w = LibXCWrapper(key)
    print >> s, '   | %s' % w.name
    for line in w.refs.split('\n'):
        print >> s, '   | *%s*' % line
    print >> s

write_if_changed('ref_functionals.rst', s.getvalue())
