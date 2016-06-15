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


import os


__all__ = ['write_rst_table', 'write_if_changed']


def write_rst_table(f, table, nhead=1):
    """Write an RST table to file f

       **Arguments:**

       f
            A writable file-like object

       table
            A list of rows. Each row must be a list of cells containing strings.

       **Optional arguments:**

       nhead
            The number of header rows
    """

    def format_cell(cell):
        if cell is None or len(cell.strip()) == 0:
            return '\ '
        else:
            return str(cell)

    # Determine the width of each column
    widths = {}
    for row in table:
        for icell, cell in enumerate(row):
            widths[icell] = max(widths.get(icell, 2), len(format_cell(cell)))

    def format_row(cells, margin):
        return ' '.join(margin + format_cell(cell).rjust(widths[icell]) + margin
                        for icell, cell in enumerate(cells))

    # construct the column markers
    markers = format_row(['='*widths[icell] for icell in xrange(len(widths))], '=')

    # top markers
    print >> f, markers

    # heading rows (if any)
    for irow in xrange(nhead):
        print >> f, format_row(table[irow], ' ')
    if nhead > 0:
        print >> f, markers

    # table body
    for irow in xrange(nhead, len(table)):
        print >> f, format_row(table[irow], ' ')

    # bottom markers
    print >> f, markers


def write_if_changed(fn, s_new):
    '''Write the string s_new to the file only if the content changes

       This is just a helper function that only really writes files when needed.
       The purpose is to avoid that Sphinx rebuilds pages that have not really
       changed.

       **Arguments:**

       fn
            The filename

       s_new
            The contents to be written to the file.
    '''
    if os.path.isfile(fn):
        # read the entire file
        with open(fn) as f:
            s_old = f.read()
        if s_new == s_old:
            print 'File %s needs no update. Skipping.' % fn
            return

    # write the new file to dis
    print 'Writing new or updated %s' % fn
    with open(fn, 'w') as f:
        f.write(s_new)
