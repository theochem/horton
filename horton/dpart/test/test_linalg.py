# -*- coding: utf-8 -*-
# Horton is a Density Functional Theory program.
# Copyright (C) 2011-2012 Toon Verstraelen <Toon.Verstraelen@UGent.be>
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


import numpy as np
from cPickle import loads

from horton import *


def test_case1():
    A = loads('cnumpy.core.multiarray\n_reconstruct\np1\n(cnumpy\nndarray\np2\n(I0\ntS\'b\'\ntRp3\n(I1\n(I6\nI6\ntcnumpy\ndtype\np4\n(S\'f8\'\nI0\nI1\ntRp5\n(I3\nS\'<\'\nNNNI-1\nI-1\nI0\ntbI00\nS\'\\xf7\\xf3\\xba4\\xf7sR@\\x89\\xbduvu\\xf4R@`\\x8f\\r\\x11\\xb5IS@{\\xe1P\\x9c\\xdb|S@\\xaabWJ\\xbc\\x90S@g\\x1bNh\\xf4\\x9aS@\\x89\\xbduvu\\xf4R@&\\xed\\x14\\xc9\\xce`S@\\xd4\\xdf\\xd7\\xa3\\xc8\\xa8S@\\xed\\x0b\\xfc\\xca\\x99\\xd4S@\\xac\\xbar5\\r\\xe7S@U\\x7f\\x0b\\x9f\\xcf\\xf1S@`\\x8f\\r\\x11\\xb5IS@\\xd4\\xdf\\xd7\\xa3\\xc8\\xa8S@\\xd9\\xd5\\xf1\\x7fn\\xe8S@\\x12\\x01\\xd9O\\x13\\x10T@\\xf0`\\xc2\\xd2C"T@F\\xe4\\t_\\xd4-T@{\\xe1P\\x9c\\xdb|S@\\xed\\x0b\\xfc\\xca\\x99\\xd4S@\\x12\\x01\\xd9O\\x13\\x10T@\\xed\\xa94\\x91\\x196T@\\xe5\\xd5\\xfb\\xa7\\xdfHT@\\xa0\\xd9\\xda\\x8bdUT@\\xaabWJ\\xbc\\x90S@\\xac\\xbar5\\r\\xe7S@\\xf0`\\xc2\\xd2C"T@\\xe5\\xd5\\xfb\\xa7\\xdfHT@\\x173\\xcd\\xa0\\xc1\\\\T@8^\\xda\\xbf5jT@g\\x1bNh\\xf4\\x9aS@U\\x7f\\x0b\\x9f\\xcf\\xf1S@F\\xe4\\t_\\xd4-T@\\xa0\\xd9\\xda\\x8bdUT@8^\\xda\\xbf5jT@\\x05r\\x1a\\x9fexT@\'\ntb.')
    b = loads("cnumpy.core.multiarray\n_reconstruct\np1\n(cnumpy\nndarray\np2\n(I0\ntS'b'\ntRp3\n(I1\n(I6\ntcnumpy\ndtype\np4\n(S'f8'\nI0\nI1\ntRp5\n(I3\nS'<'\nNNNI-1\nI-1\nI0\ntbI00\nS'\\xfb\\xa9\\x06{\\xab\\xaaS@\\xe3k\\xecv\\xea\\xffS@\\xd4\\x98f\\x87\\xe0:T@,[\\x16\\xf9\\xcfaT@&\\xa6\\x08/YvT@I\\xf3\\x1e\\xf0J\\x84T@'\ntb.")
    binding_eq = loads("(lp1\n(cnumpy.core.multiarray\n_reconstruct\np2\n(cnumpy\nndarray\np3\n(I0\ntS'b'\ntRp4\n(I1\n(I6\ntcnumpy\ndtype\np5\n(S'i8'\nI0\nI1\ntRp6\n(I3\nS'<'\nNNNI-1\nI-1\nI0\ntbI00\nS'\\x04\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x05\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x06\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x07\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x08\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\t\\x00\\x00\\x00\\x00\\x00\\x00\\x00'\ntbcnumpy.core.multiarray\nscalar\np7\n(g5\n(S'f8'\nI0\nI1\ntRp8\n(I3\nS'<'\nNNNI-1\nI-1\nI0\ntbS'G\\xf2\\xc2\\xfb\\x04\\xd6 @'\ntRp9\ntp10\na.")
    binding_ineq = loads("(lp1\n(cnumpy.core.multiarray\n_reconstruct\np2\n(cnumpy\nndarray\np3\n(I0\ntS'b'\ntRp4\n(I1\n(I6\ntcnumpy\ndtype\np5\n(S'f8'\nI0\nI1\ntRp6\n(I3\nS'<'\nNNNI-1\nI-1\nI0\ntbI00\nS'\\x00\\x00\\x00\\x00\\x00\\x00\\xf0?\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00'\ntbI0\ntp7\na(g2\n(g3\n(I0\ntS'b'\ntRp8\n(I1\n(I6\ntg6\nI00\nS'\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\xf0?\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00'\ntbI0\ntp9\na(g2\n(g3\n(I0\ntS'b'\ntRp10\n(I1\n(I6\ntg6\nI00\nS'\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\xf0?\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00'\ntbI0\ntp11\na(g2\n(g3\n(I0\ntS'b'\ntRp12\n(I1\n(I6\ntg6\nI00\nS'\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\xf0?\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00'\ntbI0\ntp13\na(g2\n(g3\n(I0\ntS'b'\ntRp14\n(I1\n(I6\ntg6\nI00\nS'\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\xf0?\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00'\ntbI0\ntp15\na(g2\n(g3\n(I0\ntS'b'\ntRp16\n(I1\n(I6\ntg6\nI00\nS'\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\xf0?'\ntbI0\ntp17\na.")

    try:
        quadratic_solver(A, b, binding_eq, binding_ineq, rcond=0)    
        assert False
    except ValueError:
        pass
