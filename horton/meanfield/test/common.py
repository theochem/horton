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


from horton import *


__all__ = ['get_some_grid', 'check_cubic_cs_wrapper', 'check_cubic_os_wrapper']


def get_some_grid(sys):
    rtf = ExpRTransform(2e-4, 2e1, 100)
    rgrid = RadialGrid(rtf)
    return BeckeMolGrid(sys, (rgrid, 110), random_rotate=False)


def check_cubic_cs_wrapper(ham, dm0, dm1, do_plot=False):
    wfn = ham.system.wfn
    fock = ham.system.lf.create_one_body()

    # evaluate stuff at dm0
    ham.invalidate()
    wfn.invalidate()
    wfn.update_dm('alpha', dm0)
    e0 = ham.compute()
    fock.reset()
    ham.compute_fock(fock, None)
    ev_00 = fock.expectation_value(dm0)
    ev_01 = fock.expectation_value(dm1)
    g0 = 2*(ev_01 - ev_00)

    # evaluate stuff at dm1
    ham.invalidate()
    wfn.invalidate()
    wfn.update_dm('alpha', dm1)
    e1 = ham.compute()
    fock.reset()
    ham.compute_fock(fock, None)
    ev_10 = fock.expectation_value(dm0)
    ev_11 = fock.expectation_value(dm1)
    g1 = 2*(ev_11 - ev_10)

    check_cubic_cs(ham, dm0, dm1, e0, e1, g0, g1, do_plot)


def check_cubic_os_wrapper(ham, dma0, dmb0, dma1, dmb1, do_plot=False):
    wfn = ham.system.wfn
    focka = ham.system.lf.create_one_body()
    fockb = ham.system.lf.create_one_body()

    # evaluate stuff at 0
    ham.invalidate()
    wfn.invalidate()
    wfn.update_dm('alpha', dma0)
    wfn.update_dm('beta', dmb0)
    e0 = ham.compute()
    focka.reset()
    fockb.reset()
    ham.compute_fock(focka, fockb)
    ev_00 = focka.expectation_value(dma0) + fockb.expectation_value(dmb0)
    ev_01 = focka.expectation_value(dma1) + fockb.expectation_value(dmb1)
    g0 = (ev_01 - ev_00)

    # evaluate stuff at 1
    ham.invalidate()
    wfn.invalidate()
    wfn.update_dm('alpha', dma1)
    wfn.update_dm('beta', dmb1)
    e1 = ham.compute()
    focka.reset()
    fockb.reset()
    ham.compute_fock(focka, fockb)
    ev_10 = focka.expectation_value(dma0) + fockb.expectation_value(dmb0)
    ev_11 = focka.expectation_value(dma1) + fockb.expectation_value(dmb1)
    g1 = (ev_11 - ev_10)

    check_cubic_os(ham, dma0, dmb0, dma1, dmb1, e0, e1, g0, g1, do_plot)
