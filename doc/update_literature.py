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


from cStringIO import StringIO

from horton.log import Biblio
from horton.context import context

from common import write_if_changed


biblio = Biblio(context.get_fn('references.bib'))

def key(item):
    return int(item[1].tags['year']), item[0]

f = StringIO()
print >> f, 'Literature'
print >> f, '##########'
items = biblio._records.items()
items.sort(key=key)
for key, reference in items:
    print >> f
    print >> f, '.. [%s] %s' % (key, reference.format_rst())
s = f.getvalue()
write_if_changed('tech_ref_literature.rst', s)
