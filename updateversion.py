#!/usr/bin/env python

import re, sys



rules = [
    ('setup.py', '^    version=\'(...+)\',$'),
    ('horton/__init__.py', '^__version__ = \'(...+)\'$'),
    ('horton/log.py', '^.* Welcome to Horton (...+)!$'),
    ('doc/conf.py', '^version = \'(...+)\'$'),
    ('doc/conf.py', '^release = \'(...+)\'$'),
    ('doc/tut_getting_started.rst', '^    http://users.ugent.be/~tovrstra/horton/horton-(...+).tar.gz.$'),
    ('doc/tut_getting_started.rst', '^    wget http://users.ugent.be/~tovrstra/horton/horton-(...+).tar.gz$'),
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
        for i in xrange(len(lines)):
            line = lines[i]
            m = r.match(line)
            if m is not None:
                lines[i] = line[:m.start(1)] + newversion + line[m.end(1):]
        with open(fn, 'w') as f:
            f.writelines(lines)
