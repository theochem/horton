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
'''Tools for avoiding recomputation and reallocation of earlier results

   In principle, the ``JustOnceClass`` and the ``Cache`` can be used
   independently, but in some cases it makes a lot of sense to combine them.
   See for example the density partitioning code in ``horton.dpart``.
'''


import numpy as np, h5py as h5, os
from horton.log import log


__all__ = ['JustOnceClass', 'just_once', 'Cache', 'Scratch']


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



class CacheItem(object):
    def __init__(self, value, own=False):
        self._value = value
        self._valid = True
        self._own = own

    @classmethod
    def from_alloc(cls, alloc):
        from horton.matrix import LinalgFactory
        if not hasattr(alloc, '__len__'):
            alloc = (alloc,)
        if isinstance(alloc[0], LinalgFactory):
            if len(alloc) < 2:
                raise TypeError('Add least one extra parameter needed to initialize a linalg thing')
            if alloc[1] == 'one_body':
                return cls(alloc[0].create_one_body(*alloc[2:]))
            elif alloc[1] == 'expansion':
                return cls(alloc[0].create_expansion(*alloc[2:]))
            else:
                raise TypeError('Not supported: %s.' % alloc[1])
        else:
            # initialize a floating point array
            array = np.zeros(alloc, float)
            log.mem.announce(array.nbytes)
            return cls(array, own=True)

    def __del__(self):
        if self._own and log is not None:
            assert isinstance(self._value, np.ndarray)
            log.mem.denounce(self._value.nbytes)

    def check_alloc(self, alloc):
        from horton.matrix import LinalgFactory
        if not hasattr(alloc, '__len__'):
            alloc = (alloc,)
        if isinstance(alloc[0], LinalgFactory):
            if len(alloc) < 2:
                raise TypeError('Add least one extra parameter needed to initialize a linalg thing')
            if alloc[1] == 'one_body':
                if not (len(alloc) == 2 or (len(alloc) == 3 and alloc[2] == self._value.nbasis)):
                    raise TypeError('The requested one-body operator is not compatible with the cached one.')
            elif alloc[1] == 'expansion':
                if not (len(alloc) == 2 or (len(alloc) == 3 and alloc[2] == self._value.nbasis) or
                        (len(alloc) == 4 and alloc[2] == self._value.nbasis and alloc[3] == self._value.nfn)):
                    raise TypeError('The requested expansion is not compatible with the cached one.')
            else:
                raise TypeError('For the moment, only one_body stuff is supported.')
        else:
            # assume a floating point array
            if not (isinstance(self._value, np.ndarray) and
                    self._value.shape == tuple(alloc) and
                    issubclass(self._value.dtype.type, float)):
                raise TypeError('The stored item does not match the given alloc.')

    def _get_value(self):
        if not self._valid:
            raise ValueError('This cached item is not valid.')
        return self._value

    value = property(_get_value)

    def _get_valid(self):
        return self._valid

    valid = property(_get_valid)

    def _get_resettable(self):
        from horton.matrix import OneBody, Expansion
        return isinstance(self._value, np.ndarray) or \
               isinstance(self._value, OneBody) or \
               isinstance(self._value, Expansion)

    resettable = property(_get_resettable)

    def invalidate(self):
        self._valid = False
        self.reset()

    def reset(self):
        from horton.matrix import OneBody
        if isinstance(self._value, np.ndarray):
            self._value[:] = 0.0
        elif isinstance(self._value, OneBody):
            self._value.reset()
        else:
            raise TypeError('Do not know how to reset %s.' % self._value)


class NoDefault(object):
    pass


no_default = NoDefault()


class Cache(object):
    '''Object that stores previously computed results.

       This is geared towards storing numerical results. The load method
       may be used to initialize new floating point arrays if they are not
       cached yet.
    '''
    def __init__(self):
        self._store = {}

    def invalidate(self):
        for key in self._store.keys():
            item = self._store[key]
            if item.resettable:
                # avoid re-allocation
                item.invalidate()
            else:
                del self._store[key]

    def load(self, *key, **kwargs):
        # check key
        if len(key) == 0:
            raise TypeError('At least one non-keyword argument is required')

        # parse kwargs
        if len(kwargs) == 0:
            alloc = None
            default = no_default
        elif len(kwargs) == 1:
            name, value = kwargs.items()[0]
            if name == 'alloc':
                alloc = value
                default = no_default
            elif name == 'default':
                alloc = None
                default = value
            else:
                raise TypeError('Only one keyword argument is allowed: alloc or default')
        else:
            raise TypeError('Only one keyword argument is allowed: alloc or default')

        item = self._store.get(key)
        if item is None:
            if alloc is None and default is no_default:
                raise KeyError(key)
            elif default is no_default:
                item = CacheItem.from_alloc(alloc)
                self._store[key] = item
                return item.value, True
            else:
                return default
        if alloc is None:
            if item.valid:
                return item.value
            else:
                raise KeyError(key)
        else:
            item.check_alloc(alloc)
            new = not item.valid
            item._valid = True # as if it is newly allocated
            return item.value, new

    def has(self, *key):
        # TODO: use pythonic "in" instead
        item = self._store.get(key)
        if item is None:
            return False
        else:
            return item.valid

    def dump(self, *args, **kwargs):
        '''Store an object in the cache.

           **Arguments:**

           key1 [, key2, ...]
                The positional arguments (except for the last) are used as a key
                for the object.

           value
                The object to be stored.

           **Optional argument:**

           own
                When set to True, the cache will take care of denouncing the
                memory usage due to this array.
        '''
        own = kwargs.pop('own', False)
        if len(kwargs) > 0:
            raise TypeError('Unknown optional arguments: %s' % kwargs.keys())
        if len(args) < 2:
            raise TypeError('At least two arguments are required: key1 and value.')
        key = args[:-1]
        value = args[-1]
        item = CacheItem(value, own)
        self._store[key] = item

    def discard(self, *key):
        '''Genuinly discards a cached item from the store, no sneaky reset.'''
        if key in self._store:
            del self._store[key]

    def rename(self, old, new):
        '''Change the name of a cached object.'''
        self._store[new] = self._store[old]
        del self._store[old]


# TODO: merge Scratch with the Cache class at some point.
class Scratch(object):
    '''A disk-based cache to reduce memory usage'''
    def __init__(self, fn, fake=False):
        self._open = True
        self._fake = fake
        if isinstance(fn, basestring):
            if fake:
                self._f = h5.File(fn, driver='core', backing_store=False, libver='latest')
            else:
                self._f = h5.File(fn, 'w', libver='latest')
            self._own = True
        elif isinstance(fn, h5.File):
            self._f = fn
            self._own = False

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()

    def __del__(self):
        self.close()

    def close(self):
        if self._open:
            self._open = False
            if self._own:
                # clean up file
                fn = self._f.filename
                self._f.close()
                if not self._fake:
                    os.remove(fn)
            else:
                # clean up everything in file
                for name in self._f.keys():
                    del self._f[name]

    def _key_to_name(self, key):
        return '_'.join(str(item) for item in key)

    def load(self, *args, **kwargs):
        if len(args) < 1:
            raise TypeError('Expecting at least one argument')
        name = self._key_to_name(args)
        output = kwargs.pop('output', None)
        if len(kwargs) > 0:
            raise TypeError('Unexpected keyword argument: %s' % kwargs.popitem()[0])
        if output is None:
            return self._f[name][:]
        else:
            self._f[name].read_direct(output)

    def dump(self, *args, **kwargs):
        data = kwargs.pop('data', None)
        if len(kwargs) > 0:
            raise TypeError('Unexpected keyword argument: %s' % kwargs.popitem()[0])
        if data is None:
            if len(args) < 2:
                raise TypeError('Expecting at least two arguments')
            name = self._key_to_name(args[:-1])
            data = args[-1]
        else:
            if len(args) < 1:
                raise TypeError('Expecting at least one argument')
            name = self._key_to_name(args)
        if name in self._f:
            self._f[name].write_direct(data)
        else:
            self._f[name] = data

    def __contains__(self, key):
        name = self._key_to_name(key)
        return name in self._f

    def rename(self, key1, key2):
        name1 = self._key_to_name(key1)
        name2 = self._key_to_name(key2)
        self._f[name2] = self._f[name1]
        del self._f[name1]
