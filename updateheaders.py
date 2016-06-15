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


from glob import glob
from fnmatch import fnmatch
import os


def strip_header(lines, closing):
    # search for the header closing line, i.e. '# --\n'
    counter = 0
    found = False
    for line in lines:
        counter += 1
        if line == closing:
            found = True
            break
    if found:
        del lines[:counter]
        # If the header closing is not found, we assume it is not present.
    # if a shebang is still present, also remove it
    if len(lines) > 0 and lines[0].startswith('#!'):
        del lines[0]
    # add a header closing line
    lines.insert(0, closing)


def fix_python(fn, lines, header_lines):
    # Do not update header of file taken from other project.
    if fn.endswith('/cpplint.py'):
        return
    # check if a shebang is present
    do_shebang = lines[0].startswith('#!')
    # remove the current header
    strip_header(lines, '# --\n')
    # add new header (insert must be in reverse order)
    for hline in header_lines[::-1]:
        lines.insert(0, ('# '+hline).strip() + '\n')
    # add a source code encoding line
    lines.insert(0, '# -*- coding: utf-8 -*-\n')
    if do_shebang:
        lines.insert(0, '#!/usr/bin/env python\n')


def fix_c(fn, lines, header_lines):
    # check for an exception line
    for line in lines:
        if 'no_update_headers' in line:
            return
    # remove the current header
    strip_header(lines, '//--\n')
    # add new header (insert must be in reverse order)
    for hline in header_lines[::-1]:
        lines.insert(0, ('// '+hline).strip() + '\n')


def fix_rst(fn, lines, header_lines):
    # check for an exception line
    for line in lines:
        if 'no_update_headers' in line:
            return
    # remove the current header
    strip_header(lines, '    : --\n')
    # add an empty line after header if needed
    if len(lines[1].strip()) > 0:
      lines.insert(1, '\n')
    # add new header (insert must be in reverse order)
    for hline in header_lines[::-1]:
        lines.insert(0, ('    : '+hline).rstrip() + '\n')
    # add comment instruction
    lines.insert(0, '..\n')


def iter_subdirs(root):
    for dn, subdns, fns in os.walk(root):
        yield dn


def main():
    source_dirs = ['.', 'doc', 'data/grids', 'scripts', 'tools', 'tools/qa'] + \
        list(iter_subdirs('horton'))

    fixers = [
        ('*.py', fix_python),
        ('*.pxd', fix_python),
        ('*.pyx', fix_python),
        ('*.txt', fix_python),
        ('*.c', fix_c),
        ('*.cpp', fix_c),
        ('*.h', fix_c),
        ('*.rst', fix_rst),
        ('*.rst.template', fix_rst),
    ]

    f = open('HEADER')
    header_lines = f.readlines()
    f.close()

    for sdir in source_dirs:
        print 'Scanning', sdir
        for fn in glob(sdir + '/*.*'):
            if not os.path.isfile(fn):
                continue
            for pattern, fixer in fixers:
                if fnmatch(fn, pattern):
                    print 'Fixing  ', fn
                    with open(fn) as f:
                        lines = f.readlines()
                    fixer(fn, lines, header_lines)
                    with open(fn, 'w') as f:
                        f.writelines(lines)
                    break


if __name__ == '__main__':
    main()
