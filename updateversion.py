#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Horton is a development platform for electronic structure methods.
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
#!/usr/bin/env python

import re, sys



rules = [
    ('setup.py', '^    version=\'(...+)\',$'),
    ('horton/__init__.py', '^__version__ = \'(...+)\'$'),
    ('horton/log.py', '^.* Welcome to Horton (...+)!$'),
    ('doc/conf.py', '^version = \'(...+)\'$'),
    ('doc/conf.py', '^release = \'(...+)\'$'),
    ('doc/download_and_install_linux.rst', '^    https://github.com/theochem/horton/releases/download/(...+)/horton-(...+).tar.gz$'),
    ('doc/download_and_install_linux.rst', '^    curl -O https://github.com/theochem/horton/releases/download/(...+)/horton-(...+).tar.gz$'),
    ('doc/download_and_install_linux.rst', '^    tar -xvzf horton-(...+).tar.gz$'),
    ('doc/download_and_install_linux.rst', '^    cd horton-(...+)$'),
    ('doc/download_and_install_mac.rst', '^    https://github.com/theochem/horton/releases/download/(...+)/horton-(...+).tar.gz$'),
    ('doc/download_and_install_mac.rst', '^    curl -O https://github.com/theochem/horton/releases/download/(...+)/horton-(...+).tar.gz$'),
    ('doc/download_and_install_mac.rst', '^    tar -xvzf horton-(...+).tar.gz$'),
    ('doc/download_and_install_mac.rst', '^    cd horton-(...+)$'),
    ('doc/index.rst', '^    Horton (...+), http://theochem.github.com/horton/,$'),
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
