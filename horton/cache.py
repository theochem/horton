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
'''Avoid recomputation of earlier results and reallocation of existing arrays

   In principle, the ``JustOnceClass`` and the ``Cache`` can be used
   independently, but in some cases it makes a lot of sense to combine them.
   See for example the density partitioning code in ``horton.part``.
'''


import numpy as np, types
from horton.log import log


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

       When all results are outdated, one can call the ``clear`` method
       to forget which methods were called already.
    '''
    def __init__(self):
        self._done_just_once = set([])

    def __clear__(self):
        self.clear()

    def clear(self):
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


def _normalize_alloc(alloc):
    '''Normalize the alloc argument of the from_alloc and check_alloc methods'''
    if not hasattr(alloc, '__len__'):
        alloc = (alloc,)
    if len(alloc) == 0:
        raise TypeError('Alloc can not be an empty list')
    return alloc


def _normalize_tags(tags):
    '''Normalize the tags argument of the CacheItem constructor'''
    if tags is None:
        return set([])
    else:
        return set(tags)


class CacheItem(object):
    '''A container for an object stored in a Cache instance'''
    def __init__(self, value, own=False, tags=None):
        '''
           **Arguments:**

           value
                The object stored in this container

           **Optional arguments:**

           own
                If True, this container will denounce the memory allocated for
                the contained object. This can only be True for a numpy array.
        '''
        self._value = value
        self._valid = True
        self._own = own
        self._tags = _normalize_tags(tags)

    @classmethod
    def from_alloc(cls, alloc, tags):
        alloc = _normalize_alloc(alloc)
        if all(isinstance(i, int) for i in alloc):
            # initialize a floating point array
            array = np.zeros(alloc, float)
            log.mem.announce(array.nbytes)
            return cls(array, own=True, tags=tags)
        else:
            # initialize a new object
            return cls(alloc[0](*alloc[1:]), tags=tags)

    def __del__(self):
        if self._own and log is not None:
            assert isinstance(self._value, np.ndarray)
            log.mem.denounce(self._value.nbytes)

    def check_alloc(self, alloc):
        alloc = _normalize_alloc(alloc)
        if all(isinstance(i, int) for i in alloc):
            # check if the array has the correct shape and dtype
            if not (isinstance(self._value, np.ndarray) and
                    self._value.shape == tuple(alloc) and
                    issubclass(self._value.dtype.type, float)):
                raise TypeError('The stored item does not match the given alloc.')
        else:
            # check if the object was initialized with compatible arguments
            try:
                if isinstance(alloc[0], type):
                    # first argument is a class
                    alloc[0].__check_init_args__(self._value, *alloc[1:])
                elif isinstance(alloc[0], types.MethodType):
                    # first argument is something else, assuming a method of a factory class
                    factory = alloc[0].__self__
                    alloc[0].__check_init_args__(factory, self._value, *alloc[1:])
                else:
                    raise NotImplementedError
            except AssertionError:
                raise TypeError('The stored item does not match the given alloc.')

    def check_tags(self, tags):
        tags = _normalize_tags(tags)
        if tags != self._tags:
            raise ValueError('Tags do not match.')

    def _get_value(self):
        if not self._valid:
            raise ValueError('This cached item is not valid.')
        return self._value

    value = property(_get_value)

    def _get_valid(self):
        return self._valid

    valid = property(_get_valid)

    def _get_tags(self):
        return self._tags

    tags = property(_get_tags)

    def clear(self):
        '''Mark the item as invalid and clear the contents of the object.

           **Returns:** A boolean indicating that clearing was successful
        '''
        self._valid = False
        if isinstance(self._value, np.ndarray):
            self._value[:] = 0.0
        elif hasattr(self._value, '__clear__') and callable(self._value.__clear__):
            self._value.__clear__()
        else:
            return False
        return True


class NoDefault(object):
    pass


no_default = NoDefault()


def _normalize_key(key):
    '''Normalize the key argument(s) of the load and dump methods'''
    if hasattr(key, '__len__') and  len(key) == 0:
        raise TypeError('At least one argument needed to specify a key.')
    # upack the key if needed
    while len(key) == 1 and isinstance(key, tuple):
        key = key[0]
    return key


class Cache(object):
    '''Object that stores previously computed results.

       The cache behaves like a dictionary with some extra features that can be
       used to avoid recomputation or reallocation.
    '''
    def __init__(self):
        self._store = {}

    def clear(self, **kwargs):
        '''Clear all items in the cache

           **Optional arguments:**

           dealloc
                When set to True, the items are really removed from memory.

           tags
                Limit the items cleared to those who have at least one tag
                that matches one of the given tags. When this argument is used
                and it contains at least one tag, items with no tags are not
                cleared.
        '''
        # Parse kwargs. This forces the caller to use keywords in order to avoid
        # confusion.
        dealloc = kwargs.pop('dealloc', False)
        tags = kwargs.pop('tags', None)
        if len(kwargs) > 0:
            raise TypeError('Unexpected arguments: %s' % kwargs.keys())
        # actual work
        tags = _normalize_tags(tags)
        for key, item in self._store.items():
            if len(tags) == 0 or len(item.tags & tags) > 0:
                self.clear_item(key, dealloc=dealloc)

    def clear_item(self, *key, **kwargs):
        '''Clear a selected item from the cache

           **Optional arguments:**

           dealloc
                When set to True, the item is really removed from memory.
        '''
        key = _normalize_key(key)
        dealloc = kwargs.pop('dealloc', False)
        if len(kwargs) > 0:
            raise TypeError('Unexpected arguments: %s' % kwargs.keys())
        item = self._store.get(key)
        if item is None:
            return
        cleared = False
        if not dealloc:
            cleared = item.clear()
        if not cleared:
            del self._store[key]

    def load(self, *key, **kwargs):
        '''Get a value from the cache

           **Arguments:**

           key0 [key1 ...]
                All positional arguments are used as keys to identify the cached
                value.

           **Optional arguments:**

           alloc
                Parameters used to allocate a cached value if it is not present
                yet. This argument can take two forms. When integer or a
                tuple of integers is given, an array is allocated.
                Alternatively, a tuple may be given whose first element is a
                constructor, and further elements are arguments for that
                constructor.

           default
                A default value that is returned when the key does not exist in
                the cache. This default value is not stored in the cache.

           tags
                When alloc is used and a new object is thereby created or
                reused, it will get these tags. This argument is only allowed if
                the alloc argument is present. In case no new object is
                allocated, the given tags must match those already present.

           The optional argument alloc and default are both meant to handle
           situations when the key has not associated value. Hence they can not
           be both present.
        '''
        key = _normalize_key(key)

        # parse kwargs
        alloc = kwargs.pop('alloc', None)
        default = kwargs.pop('default', no_default)
        tags = kwargs.pop('tags', None)
        if not (alloc is None or default is no_default):
            raise TypeError('The optional arguments alloc and default can not be used at the same time.')
        if tags is not None and alloc is None:
            raise TypeError('The tags argument is only allowed when the alloc argument is present.')
        if len(kwargs) > 0:
            raise TypeError('Unknown optional arguments: %s' % kwargs.keys())

        # get the item from the store and decide what to do
        item = self._store.get(key)
        # there are three behaviors, depending on the keyword argumentsL
        if alloc is not None:
            # alloc is given. hence two return values: value, new
            if item is None:
                # allocate a new item and store it
                item = CacheItem.from_alloc(alloc, tags)
                self._store[key] = item
                return item.value, True
            elif not item.valid:
                try:
                    # try to reuse the same memroy
                    item.check_alloc(alloc)
                    item._valid = True # as if it is newly allocated
                    item.check_tags(tags)
                except TypeError:
                    # if reuse fails, reallocate
                    item = CacheItem.from_alloc(alloc, tags)
                    self._store[key] = item
                return item.value, True
            else:
                item.check_alloc(alloc)
                item.check_tags(tags)
                return item.value, False
        elif default is not no_default:
            # a default value is given, it is not stored
            if item is None or not item.valid:
                return default
            else:
                return item.value
        else:
            # no optional arguments are given
            if item is None or not item.valid:
                raise KeyError(key)
            else:
                return item.value

    def __contains__(self, key):
        key = _normalize_key(key)
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

           tags
                Tags to be associated with the object
        '''
        own = kwargs.pop('own', False)
        tags = kwargs.pop('tags', None)
        if len(kwargs) > 0:
            raise TypeError('Unknown optional arguments: %s' % kwargs.keys())
        if len(args) < 2:
            raise TypeError('At least two arguments are required: key1 and value.')
        key = _normalize_key(args[:-1])
        value = args[-1]
        item = CacheItem(value, own, tags)
        self._store[key] = item

    def __len__(self):
        return sum(item.valid for item in self._store.itervalues())

    def __getitem__(self, key):
        return self.load(key)

    def __setitem__(self, key, value):
        return self.dump(key, value)

    def __iter__(self):
        return self.iterkeys()

    def iterkeys(self, tags=None):
        '''Iterate over the keys of all valid items in the cache.'''
        tags = _normalize_tags(tags)
        for key, item in self._store.iteritems():
            if item.valid and (len(tags) == 0 or len(item.tags & tags) > 0):
                yield key

    def itervalues(self, tags=None):
        '''Iterate over the values of all valid items in the cache.'''
        tags = _normalize_tags(tags)
        for item in self._store.itervalues():
            if item.valid and (len(tags) == 0 or len(item.tags & tags) > 0):
                yield item.value

    def iteritems(self, tags=None):
        '''Iterate over all valid items in the cache.'''
        tags = _normalize_tags(tags)
        for key, item in self._store.iteritems():
            if item.valid and (len(tags) == 0 or len(item.tags & tags) > 0):
                yield key, item.value
