#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Horton is a Density Functional Theory program.
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


from glob import glob
import os

def strip_header(lines, closing):
    # search for the header closing line, i.e. '#--\n'
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
    # add a header closing line
    lines.insert(0, closing)


def fix_python(lines, header_lines):
    # check if a shebang is present
    do_shebang = lines[0].startswith('#!')
    # remove the current header
    strip_header(lines, '#--\n')
    # add new header (insert must be in reverse order)
    for hline in header_lines[::-1]:
        lines.insert(0, ('# '+hline).strip() + '\n')
    # add a source code encoding line
    lines.insert(0, '# -*- coding: utf-8 -*-\n')
    if do_shebang:
        lines.insert(0, '#!/usr/bin/env python\n')


def fix_c(lines, header_lines):
    # check for an exception line
    for line in lines:
        if 'no_update_headers' in line:
            return
    # remove the current header
    strip_header(lines, '//--\n')
    # add new header (insert must be in reverse order)
    for hline in header_lines[::-1]:
        lines.insert(0, ('// '+hline).strip() + '\n')


def main():
    source_dirs = [
        '.',
        'horton',
        'horton/test',
        'horton/grid',
        'horton/grid/test',
        'horton/gbasis',
        'horton/gbasis/test',
    ]

    fixers = [
        ('.py', fix_python),
        ('.pxd', fix_python),
        ('.pyx', fix_python),
        ('.c', fix_c),
        ('.cpp', fix_c),
        ('.h', fix_c),
    ]

    f = open('HEADER')
    header_lines = f.readlines()
    f.close()

    for sdir in source_dirs:
        print 'Scanning', sdir
        for fn in glob(sdir + '/*.*'):
            if not os.path.isfile(fn):
                continue
            for ext, fixer in fixers:
                if fn.endswith(ext):
                    print 'Fixing  ', fn
                    f = file(fn)
                    lines = f.readlines()
                    f.close()
                    fixer(lines, header_lines)
                    f = file(fn, 'w')
                    f.writelines(lines)
                    f.close()
                    break


if __name__ == '__main__':
    main()
