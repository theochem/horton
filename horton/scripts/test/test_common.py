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


import h5py as h5, argparse

from horton.scripts.common import *


def test_parse_pbc():
    assert (parse_pbc('111') == [1, 1, 1]).all()
    assert (parse_pbc('000') == [0, 0, 0]).all()
    assert (parse_pbc('101') == [1, 0, 1]).all()
    assert (parse_pbc('001') == [0, 0, 1]).all()


def test_store_args():
    with h5.File('horton.scripts.test.test_common.test_store_args.h5', driver='core', backing_store=False) as f:
        namespace = argparse.Namespace()
        namespace.a = 1
        namespace.b = 'a'
        namespace.c = None
        namespace.d = 2.5
        store_args(namespace, f)
        assert f.attrs['arg_a'] == 1
        assert f.attrs['arg_b'] == 'a'
        assert 'arg_c' not in f.attrs
        assert f.attrs['arg_d'] == 2.5
        assert 'pwd' in f.attrs
        assert 'cmdline' in f.attrs
        assert 'datetime' in f.attrs
