#!/usr/bin/env python

import re, sys



rules = [
    ('setup.py', '^    version=\'(...+)\',$'),
    ('horton/__init__.py', '^__version__ = \'(...+)\'$'),
    ('horton/log.py', '^.* Welcome to Horton (...+)!$'),
    ('doc/conf.py', '^version = \'(...+)\'$'),
    ('doc/conf.py', '^release = \'(...+)\'$'),
    ('doc/tut_getting_started.rst', '^    https://github.com/theochem/horton/releases/download/(...+)/horton-(...+).tar.gz$'),
    ('doc/tut_getting_started.rst', '^    wget https://github.com/theochem/horton/releases/download/(...+)/horton-(...+).tar.gz$'),
    ('doc/tut_getting_started.rst', '^    tar -xvzf horton-(...+).tar.gz$'),
    ('doc/tut_getting_started.rst', '^    cd horton-(...+)$'),
    ('doc/tut_cite.rst', '^    Horton (...+), http://theochem.github.com/horton/,$'),
]


if __name__ == '__main__':
    newversion = sys.argv[1]

    for fn, regex in rules:
        r = re.compile(regex)
        with open(fn) as f:
            lines = f.readlines()
        for iline in xrange(len(lines)):
            line = lines[iline]
            m = r.match(line)
            if m is not None:
                for igroup in xrange(m.lastindex, 0, -1):
                    line = line[:m.start(igroup)] + newversion + line[m.end(igroup):]
                lines[iline] = line
        with open(fn, 'w') as f:
            f.writelines(lines)
