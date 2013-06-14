#!/usr/bin/env python

import re, sys



rules = [
    ('setup.py', '^    version=\'(...+)\',$'),
    ('horton/__init__.py', '^__version__ = \'(...+)\'$'),
    ('doc/conf.py', '^version = \'(...+)\'$'),
    ('doc/conf.py', '^release = \'(...+)\'$'),
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
