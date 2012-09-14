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
'''Tools for avoiding recomputation of earlier results

   In principle, the ``JustOnceClass`` and the ``Cache`` can be used
   independently, but in some cases it makes a lot of sense to combine them.
   See for example. the density partitioning code in ``horton.dpart``.
'''


import numpy as np


__all__ = ['JustOnceClass', 'just_once', 'Cache']


class JustOnceClass(object):
    '''Base class for classes with methods that should never be executed twice.

       In typically applications, these methods get called many times, but
       only during the first call, an actual computation is carried out. This
       way, the caller can safely call a method, just to make sure that a
       required result is computed.

       All methods in the subclasses that should have this feature, must be
       given the ``just_once`` decoratore, e.g. ::

           class Example(JustOnceClass):
               @just_once
               def do_something():
                   self.foo = self.bar

       When all results are outdate, one can call the ``invalidate`` method
       to forget which methods were called already.
    '''
    def __init__(self):
        self._done_just_once = set([])

    def invalidate(self):
        self._done_just_once = set([])


def just_once(fn):
    def wrapper(instance):
        if not hasattr(instance, '_done_just_once'):
            raise TypeError('Missing hidden _done_just_once. Forgot to call JustOnceClass.__init__()?')
        if fn.func_name in instance._done_just_once:
            return
        fn(instance)
        instance._done_just_once.add(fn.func_name)
    wrapper.__doc__ = fn.__doc__
    return wrapper


class Cache(object):
    '''Object that stores previously computed results.

       This is geared towards storing numerical results. The load method
       may be used to initialize new floating point arrays if they are not
       cached yet.
    '''
    def __init__(self):
        self._store = {}

    def invalidate():
        self._store = {}

    def load(self, *key, **kwargs):
        # check key
        if len(key) == 0:
            raise TypeError('At least one non-keyword argument is required')

        # parse kwargs
        if len(kwargs) == 0:
            newshape = None
        elif len(kwargs) == 1:
            name, newshape = kwargs.items()[0]
            if name != 'newshape':
                raise TypeError('Only one keyword argument is allowed: newshape')
        else:
            raise TypeError('Only one keyword argument is allowed: newshape')

        value = self._store.get(key)
        if value is None:
            if newshape is None:
                raise KeyError('Could not find item %s' % repr(key))
            else:
                value = np.zeros(newshape, float)
                self._store[key] = value
                return value, True
        if newshape is None:
            return value
        else:
            if not hasattr(newshape, '__len__'):
                newshape = (newshape,)
            if not (isinstance(value, np.ndarray) and value.shape == newshape and issubclass(value.dtype.type, float)):
                raise TypeError('The stored result does not have the excpected type (newshape).')
            return value, False

    def dump(self, *args):
        if len(args) < 2:
            raise TypeError('At least two arguments are required: key1 and value.')
        key = args[:-1]
        value = args[-1]
        self._store[key] = value
