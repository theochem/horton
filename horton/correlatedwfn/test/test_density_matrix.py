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


import numpy as np
from nose.tools import assert_raises
from horton import *

geminal = np.array([[ 1.92138675e-06,-4.28103597e-04,-4.28103597e-04,-1.53442827e-03],
                    [ 8.85797056e-05,-4.50696243e-04,-4.50696243e-04,-1.58078703e-03],
                    [-1.58271023e-02,-8.62977027e-03,-8.62977027e-03,-2.74604281e-02],
                    [-5.13051053e-01,-4.81018842e-02,-4.81018842e-02, 2.41968088e-03],
                    [-5.41777471e-03,-9.02405302e-03,-1.91841957e-01,-9.36965909e-03],
                    [-5.41777471e-03,-1.91841957e-01,-9.02405302e-03,-9.36965909e-03]])
lagrange = np.array([[-2.26329664e-05,-4.41712415e-04,-4.41712415e-04,-1.53417147e-03],
                     [ 4.70936832e-05,-4.62497383e-04,-4.62497383e-04,-1.58047677e-03],
                     [-1.43637183e-02,-9.25111933e-03,-9.25111933e-03,-2.75683948e-02],
                     [-4.07653846e-01,-3.26007215e-02,-3.26007215e-02, 1.74851520e-03],
                     [ 1.00308707e-02,-9.35939187e-03,-1.85828494e-01,-9.51688460e-03],
                     [ 1.00308707e-02,-1.85828494e-01,-9.35939187e-03,-9.51688460e-03]])

def test_response_3dm_uuu():
    lf = DenseLinalgFactory()
    onebody1 = lf.create_one_body(6,4)
    onebody2 = lf.create_one_body(6,4)
    threebody = lf.create_three_body(10)
    onebody1.assign_array(geminal)
    onebody2.assign_array(lagrange)
    threebody.compute_response_three_dm_ap1rog(onebody1, onebody2, 'uuu')

    # c_ab = c_ai.c_ib:
    passc = False
    check = np.zeros((10,10,10))
    product = onebody1.contract_onebody(onebody2)
    for i in range(6):
        for j in range(6):
            for k in range(6):
                val = 0.0
                for a in range(4):
                    for l in range(6):
                        if l is not i and l is not k and l is not j:
                            if i is not j and i is not k:
                                if k is not j:
                                    val += geminal[l,a]*lagrange[l,a]
                                    passc = True
                if passc:
                    check[i,j,k] = 1+val-product
                    passc = False
            for a in range(4):
                val = 0.0
                for l in range(6):
                    if l is not i and l is not j:
                        if i is not j:
                            val += geminal[l,a]*lagrange[l,a]
                            passc = True
                if passc:
                    check[(a+6),i,j] = val
                    check[i,(a+6),j] = val
                    check[i,j,(a+6)] = val
                    passc = False
    assert np.allclose(check, threebody._array)

def test_response_3dm_uud():

    lf = DenseLinalgFactory()
    onebody1 = lf.create_one_body(6,4)
    onebody2 = lf.create_one_body(6,4)
    threebody = lf.create_three_body(10)
    onebody1.assign_array(geminal)
    onebody2.assign_array(lagrange)
    threebody.compute_response_three_dm_ap1rog(onebody1, onebody2, 'uud')

    # c_ab = c_ai.c_ib:
    check = np.zeros((10,10,10))
    product = onebody1.contract_onebody(onebody2)
    for i in range(6):
        for j in range(6):
            for k in range(6):
                val = 0.0
                for a in range(4):
                    for l in range(6):
                        if l is not i and l is not k and l is not j:
                            val += geminal[l,a]*lagrange[l,a]
                check[i,j,k] = 1+val-product
                if i is j:
                    check[i,i,k] = 0.0
            for a in range(4):
                val = 0.0
                for l in range(6):
                    if l is not i and l is not j:
                        val += geminal[l,a]*lagrange[l,a]
                check[(a+6),i,j] = val
                check[i,(a+6),j] = val
                if i is not j:
                    check[i,j,(a+6)] = val
        for a in range(4):
            val = 0.0
            for l in range(6):
                if l is not i:
                    val += geminal[l,a]*lagrange[l,a]
            check[(a+6),i,(a+6)] = val
            check[i,(a+6),(a+6)] = val
    assert np.allclose(check, threebody._array)

def test_response_3dm_uud_offdiagonal():

    lf = DenseLinalgFactory()
    onebody1 = lf.create_one_body(6,4)
    onebody2 = lf.create_one_body(6,4)
    threebody = lf.create_three_body(10)
    onebody1.assign_array(geminal)
    onebody2.assign_array(lagrange)
    threebody.compute_response_three_dm_ap1rog(onebody1, onebody2, 'uudoff')

    # c_ab = c_ai.c_ib:
    check = np.zeros((10,10,10))
    product = onebody1.contract_onebody(onebody2)
    for i in range(6):
        for j in range(6):
            for k in range(6):
                val = 0.0
                for a in range(4):
                    val += geminal[j,a]*lagrange[k,a]
                check[i,j,k] = val
                if i is j:
                    check[i,i,k] = 0.0
                if k is i:
                    check[i,j,i] = 0.0
                if k is j:
                    check[i,j,j] = 0.0
            for a in range(4):
                check[(a+6),j,i] = geminal[j,a]*lagrange[i,a]
                check[i,(a+6),j] = lagrange[j,a]
                val = 0.0
                for l in range(6):
                    if l is not j:
                        for c in range(4):
                            if c is not a:
                                val += lagrange[l,c]*(geminal[l,c]*geminal[j,a]+geminal[l,a]*geminal[j,c])
                check[i,j,(a+6)] = geminal[j,a]*(1-product)+val
                if i is j:
                    check[i,(a+6),j] = 0.0
                    check[(a+6),j,i] = 0.0
                    check[i,j,(a+6)] = 0.0
        for a in range(4):
            for b in range(4):
                val = 0.0
                for l in range(6):
                    if l is not i:
                        val += geminal[l,b]*lagrange[l,a]
                check[i,(a+6),(b+6)] = val
        check[i,i,i] = 0.0
    assert np.allclose(check, threebody._array)

def test_response_3dm_udd():

    lf = DenseLinalgFactory()
    onebody1 = lf.create_one_body(6,4)
    onebody2 = lf.create_one_body(6,4)
    threebody = lf.create_three_body(10)
    onebody1.assign_array(geminal)
    onebody2.assign_array(lagrange)
    threebody.compute_response_three_dm_ap1rog(onebody1, onebody2, 'udd')

    # c_ab = c_ai.c_ib:
    check = np.zeros((10,10,10))
    product = onebody1.contract_onebody(onebody2)
    for i in range(6):
        for j in range(6):
            for k in range(6):
                val = 0.0
                for a in range(4):
                    for l in range(6):
                        if l is not i and l is not k and l is not j:
                            val += geminal[l,a]*lagrange[l,a]
                check[i,j,k] = 1+val-product
                if k is j:
                    check[i,k,k] = 0.0
            for a in range(4):
                val = 0.0
                for l in range(6):
                    if l is not i and l is not j:
                        val += geminal[l,a]*lagrange[l,a]
                if i is not j:
                   check[(a+6),i,j] = val
            for a in range(4):
                val = 0.0
                for l in range(6):
                    if l is not j and l is not i:
                        val += geminal[l,a]*lagrange[l,a]
                check[i,(a+6),j] = val
                check[i,j,(a+6)] = val
        for a in range(4):
            val = 0.0
            for l in range(6):
                if l is not i:
                    val += geminal[l,a]*lagrange[l,a]
            check[(a+6),i,(a+6)] = val
            check[i,(a+6),(a+6)] = val
    assert np.allclose(check, threebody._array)

def test_response_3dm_udd_offdiagonal():

    lf = DenseLinalgFactory()
    onebody1 = lf.create_one_body(6,4)
    onebody2 = lf.create_one_body(6,4)
    threebody = lf.create_three_body(10)
    onebody1.assign_array(geminal)
    onebody2.assign_array(lagrange)
    threebody.compute_response_three_dm_ap1rog(onebody1, onebody2, 'uddoff')

    # c_ab = c_ai.c_ib:
    check = np.zeros((10,10,10))
    product = onebody1.contract_onebody(onebody2)
    for i in range(6):
        for j in range(6):
            for k in range(6):
                val = 0.0
                for a in range(4):
                    val += geminal[i,a]*lagrange[j,a]
                check[i,j,k] = val
                if i is j:
                    check[i,i,k] = 0.0
                if k is i:
                    check[i,j,i] = 0.0
                if k is j:
                    check[i,j,j] = 0.0
            for a in range(4):
                check[(a+6),i,j] = lagrange[i,a]
                val = 0.0
                for l in range(6):
                    if l is not i:
                        for c in range(4):
                            if c is not a:
                                val += lagrange[l,c]*(geminal[l,c]*geminal[i,a]+geminal[l,a]*geminal[i,c])
                check[i,(a+6),j] = geminal[i,a]*(1-product)+val
                if i is j:
                    check[i,(a+6),j] = 0.0
                    check[(a+6),j,i] = 0.0
                    check[i,j,(a+6)] = 0.0
        for a in range(4):
            for b in range(4):
                val = 0.0
                for l in range(6):
                    if l is not i:
                        val += geminal[l,b]*lagrange[l,a]
                check[(a+6),(b+6),i] = val
        check[i,i,i] = 0.0
    assert np.allclose(check, threebody._array)

def test_response_4dm_udud():

    lf = DenseLinalgFactory()
    onebody1 = lf.create_one_body(6,4)
    onebody2 = lf.create_one_body(6,4)
    onebody1.assign_array(geminal)
    onebody2.assign_array(lagrange)
    onedm = lf.create_zero_body(10)
    twodm = lf.create_one_body(10)
    onedm.compute_response_one_dm_ap1rog(onebody1, onebody2)
    twodm.compute_response_two_dm_ap1rog(onedm, onebody1, onebody2, 'pqpq')
    fourbody = lf.create_one_body(10)
    fourbody.compute_response_four_dm_ap1rog(twodm, onebody1, onebody2, 'udud')

    # c_ab = c_ai.c_ib:
    passc = False
    check = np.zeros((10,10))
    product = onebody1.contract_onebody(onebody2)
    for i in range(6):
        for j in range(6):
            val = 0.0
            for a in range(4):
                for l in range(6):
                    if l is not i and l is not j:
                        if i is not j:
                            val += geminal[l,a]*lagrange[l,a]
                            passc = True
            if passc:
                check[i,j] = 1+val-product
                passc = False
        for a in range(4):
            val = 0.0
            for l in range(6):
                if l is not i:
                    val += geminal[l,a]*lagrange[l,a]
                    passc = True
            if passc:
                check[(a+6),i] = val
                check[i,(a+6)] = val
                passc = False
    assert np.allclose(check, fourbody._array)
