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

# sort the functions

splitkeys = []
for key in keys:
    words = key.split('_')
    if words[0] == 'hyb':
        prefix = '_'.join(words[:2])
        mid = words[2]
        suffix = '_'.join(words[3:])
    else:
        prefix = words[0]
        mid = words[1]
        suffix = '_'.join(words[2:])
    splitkeys.append((prefix, mid, suffix))

def cmp_prefix(prefix1, prefix2):
    l = ['lda', 'gga', 'hyb_gga']
    pos1 = l.index(prefix1)
    pos2 = l.index(prefix2)
    return cmp(pos1, pos2)

def cmp_middle(middle1, middle2):
    l = ['k', 'x', 'c', 'xc']
    pos1 = l.index(middle1)
    pos2 = l.index(middle2)
    return cmp(pos1, pos2)

def cmp_splitkey(sk1, sk2):
    result = cmp_prefix(sk1[0], sk2[0])
    if result == 0:
        result = cmp_middle(sk1[1], sk2[1])
    if result == 0:
        result = cmp(sk1[2], sk2[2])
    return result

splitkeys.sort(cmp=cmp_splitkey)
keys = []
for splitkey in splitkeys:
    splitkey = [part for part in splitkey if len(part) > 0]
    keys.append('_'.join(splitkey))

# make a rst table of all functionals

from horton import RestrictedLibXCWrapper
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
    w = RestrictedLibXCWrapper(key)
    print >> s, '   | %s' % w.name
    for line in w.refs.split('\n'):
        print >> s, '   | *%s*' % line
    print >> s

write_if_changed('ref_functionals.rst', s.getvalue())
