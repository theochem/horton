# -*- coding: utf-8 -*-
# HORTON: Helpful Open-source Research TOol for N-fermion systems.
# Copyright (C) 2011-2019 The HORTON Development Team
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
'''Utility functions'''


__all__ = ['typecheck_geo', 'check_type', 'check_options', 'doc_inherit']


def typecheck_geo(coordinates=None, numbers=None, pseudo_numbers=None,
                  need_coordinates=True, need_numbers=True,
                  need_pseudo_numbers=True):
    '''Type check a molecular geometry specification

       **Arguments:**

       coordinates
            A (N, 3) float array with Cartesian coordinates of the atoms.

       numbers
            A (N,) int vector with the atomic numbers.

       **Optional arguments:**

       pseudo_numbers
            A (N,) float array with pseudo-potential core charges.

       need_coordinates
            When set to False, the coordinates can be None, are not type checked
            and not returned.

       need_numbers
            When set to False, the numbers can be None, are not type checked
            and not returned.

       need_pseudo_numbers
            When set to False, the pseudo_numbers can be None, are not type
            checked and not returned.

       **Returns:** ``[natom]`` + all arguments that were type checked. The
       pseudo_numbers argument is converted to a floating point array.
    '''
    # Determine natom
    if coordinates is not None:
        natom = len(coordinates)
    elif numbers is not None:
        natom = len(numbers)
    elif pseudo_numbers is not None:
        natom = len(pseudo_numbers)
    else:
        raise TypeError('At least one argument is required and should not be None')

    # Typecheck coordinates:
    if coordinates is None:
        if need_coordinates:
            raise TypeError('Coordinates can not be None.')
    else:
        if coordinates.shape != (natom, 3) or not issubclass(coordinates.dtype.type, float):
            raise TypeError('The argument centers must be a float array with shape (natom,3).')

    # Typecheck numbers
    if numbers is None:
        if need_numbers:
            raise TypeError('Numbers can not be None.')
    else:
        if numbers.shape != (natom,) or not issubclass(numbers.dtype.type, int):
            raise TypeError('The argument numbers must be a vector with length natom.')

    # Typecheck pseudo_numbers
    if pseudo_numbers is None:
        if need_pseudo_numbers:
            pseudo_numbers = numbers.astype(float)
    else:
        if pseudo_numbers.shape != (natom,):
            raise TypeError('The argument pseudo_numbers must be a vector with length natom.')
        if not issubclass(pseudo_numbers.dtype.type, float):
            pseudo_numbers = pseudo_numbers.astype(float)

    # Collect return values
    result = [natom, ]
    if need_coordinates:
        result.append(coordinates)
    if need_numbers:
        result.append(numbers)
    if need_pseudo_numbers:
        result.append(pseudo_numbers)
    return result


def check_type(name, instance, *Classes):
    '''Check type of argument with given name against list of types

       **Arguments:**

       name
            The name of the argument being checked.

       instance
            The object being checked.

       Classes
            A list of allowed types.
    '''
    if len(Classes) == 0:
        raise TypeError('Type checking with an empty list of classes. This is a simple bug!')
    match = False
    for Class in Classes:
        if isinstance(instance, Class):
            match = True
            break
    if not match:
        classes_parts = ['\'', Classes[0].__name__, '\'']
        for Class in Classes[1:-1]:
            classes_parts.extend([', ``', Class.__name__, '\''])
        if len(Classes) > 1:
            classes_parts.extend(['or \'', Class.__name__, '\''])
        raise TypeError('The argument \'%s\' must be an instance of %s. Got a \'%s\' instance instead.' % (
            name, ''.join(classes_parts), instance.__class__.__name__
        ))


def check_options(name, select, *options):
    '''Check if a select is in the list of options. If not raise ValueError

       **Arguments:**

       name
            The name of the argument.

       select
            The value of the argument.

       options
            A list of allowed options.
    '''
    if select not in options:
        formatted = ', '.join(['\'%s\'' % option for option in options])
        raise ValueError('The argument \'%s\' must be one of: %s' % (name, formatted))


def doc_inherit(base_class):
    """Docstring inheriting method descriptor

       doc_inherit decorator

       Usage:

       .. code-block:: python

            class Foo(object):
                def foo(self):
                    "Frobber"
                    pass

            class Bar(Foo):
                @doc_inherit(Foo)
                def foo(self):
                    pass

       Now, ``Bar.foo.__doc__ == Bar().foo.__doc__ == Foo.foo.__doc__ ==
       "Frobber"``
    """
    def decorator(method):
        overridden = getattr(base_class, method.__name__, None)
        if overridden is None:
            raise NameError('Can\'t find method \'%s\' in base class.')
        method.__doc__ = overridden.__doc__
        return method
    return decorator
