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


import os


__all__ = ['write_if_changed']


def write_if_changed(fn, s_new):
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
