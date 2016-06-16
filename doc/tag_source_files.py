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
'''Script that inserts a commit id in the source files that accompany the html files

   This is useful when people download the source files and start editing and
   e-mailing them instead of using GIT.
'''

import subprocess, os, sys
from fnmatch import fnmatch


def iter_txt_files(root):
    for dn, subdns, fns in os.walk(root):
        for fn in fns:
            if fnmatch(fn, '*.txt'):
                yield os.path.join(dn, fn)


def main(root):
    from conf import release

    print 'Tagging files in %s with release %s' % (root, release)

    for fn_txt in iter_txt_files(root):
        with open(fn_txt) as f:
            lines = f.readlines()
        if lines[1].startswith('    DOCUMENTATION BUILT FROM'):
            lines = lines[3:]
        lines.insert(0, '\n')
        lines.insert(0, '    DOCUMENTATION BUILT FROM RELEASE: %s\n' % release)
        lines.insert(0, '.. \n')
        with open(fn_txt, 'w') as f:
            f.writelines(lines)

if __name__ == '__main__':
    main(sys.argv[1])
