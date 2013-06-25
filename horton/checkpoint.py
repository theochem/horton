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
"""This field defines all the parts of the System class than can be written
   to a checkpoint file, and how that should be done in each case.

   The philosophy of this module is to take away as much of the checkpointing
   machinery out of the rest of Horton. The only (somewhat unavoidable)
   exceptions to this philosophy are the ``load_checkpoint`` function in
   ``horton.io``, the ``from_hdf5`` and ``to_hdf5`` methods in several parts
   of Horton and calls to System.update_chk in several parts of the code.

   TODO: (long future) make the fields clever enough such that they can convert
   data from one linalgfactory into another.
"""


import h5py as h5, numpy as np

from horton.cache import Cache


__all__ = ['attribute_register', 'load_hdf5_low', 'dump_hdf5_low']


class CHKField(object):
    def __init__(self, att_name, key=None, tags=None):
        """
           **Argument:**

           att_name
                The name of the attribute of this field

           **Optional argument:**

           key
                If the attribute is a dictionary, the key refers to an element
                in the dictionary. This must be a string

           tags
                These may be provided for items in a Cache object.

           The attribute of the System class may be None, a numpy array or
           an object that has to_hdf5 and from hdf5 methods method. In the last
           case, the HDF5 group must have a 'class' attribute.
        """
        if not isinstance(att_name, basestring):
            raise TypeError('att_name must be a string')
        if key is not None and not isinstance(key, basestring):
            raise TypeError('key must be a string, when given')
        self.att_name = att_name
        self.key = key
        self.tags = tags

    def read(self, chk, lf):
        """Read and return an the attribute from the chk file

           **Argument:**

           chk
                A h5.File object corresponding to a Horton checkpoint file.

           lf
                A LinalgFractry class. For now, this is just used as a double
                check.

           Depending on the contents of the checkpoint file corresponding to
           this field, the reading is done in a specific way. When not data is
           present, None is returned. When a dataset is present, a Numpy array
           is returned. When a group is found, the 'class' attribute is used
           to create an object of the right class through its ``from_hdf5``
           class method.
        """
        # A) find the right item in the checkpoint file. If it can not be found,
        #    None is returned.
        item = chk.get(self.att_name)
        if item is None:
            return None
        if self.key is not None:
            item = item.get(self.key)
            if item is None:
                return None

        # B) Construct the corresponding Python object
        return load_hdf5_low(item, lf)

    def write(self, chk, system):
        """Write an attribute to a Horton checkpoint file.

           **Arguments:**

           chk
                A h5.File object corresponding to a Horton checkpoint file.

           system
                A System instance.

           Depending on the attribute, the writing to the checkpoint file is
           done in different ways. When the attribute is None, the corresponding
           dataset of group is deleted. When it is a numpy array, it will be
           written to the checkpoint file as a dataset. When it is of another
           type, the object must have a ``to_hdf5`` file. This method will be
           given a group in the checkpoint file to write its contents too. In
           case of composed objects, the ``to_hdf5`` method may call to_hdf5
           methods of its attributes.
        """
        # A) get the object that must be written to the checkpoint file
        att = getattr(system, self.att_name, None)
        if att is None:
            return
        if self.key is not None:
            if isinstance(att, dict):
                att = att.get(self.key)
            elif isinstance(att, Cache):
                att = att.load(self.key, default=None)
            if att is None:
                return

        # B) Get the name of the object (and the group)
        grp = chk
        name = self.att_name
        if self.key is not None:
            grp = chk.require_group(name)
            name = self.key

        # C) Actual write
        dump_hdf5_low(grp, name, att)


def load_hdf5_low(item, lf):
    '''Load an object from an h5py object from a checkpoint file

       **Arguments:**

       item
            A HD5 Dataset or group.
    '''
    if isinstance(item, h5.Dataset):
        if len(item.shape) > 0:
            # convert to a numpy array
            return np.array(item)
        else:
            # convert to a scalar
            return item[()]
    elif isinstance(item, h5.Group):
        class_name = item.attrs.get('class')
        if class_name is None:
            # assuming that an entire dictionary must be read. only
            # read datasets. raise error when group is encountered.
            result = {}
            for key, subitem in item.iteritems():
                result[key] = load_hdf5_low(subitem, lf)
            return result
        else:
            # special constructor. the class is found with the imp module
            cls = __import__('horton', fromlist=[class_name]).__dict__[class_name]
            return cls.from_hdf5(item, lf)


def dump_hdf5_low(grp, name, att):
    '''Low-level routine used to dump data to an HDF5 file.

       grp
            A HDF5 group.

       name
            The name to be used for this object.

       att
            The object to be written.
    '''
    if isinstance(att, int) or isinstance(att, float) or isinstance(att, np.ndarray):
        # Simply overwrite old dataset
        if name in grp:
            del grp[name]
        grp[name] = att
    elif isinstance(att, dict):
        # Simply overwrite old datagroup
        if name in grp:
            del grp[name]
        grp = grp.create_group(name)
        for key, value in att.iteritems():
            dump_hdf5_low(grp, key, value)
    else:
        grp = grp.require_group(name)
        # clear the group if anything was present
        for key in grp.keys():
            del grp[key]
        for key in grp.attrs.keys():
            del grp.attrs[key]
        att.to_hdf5(grp)
        # The following is needed to create object of the right type when
        # reading from the checkpoint:
        grp.attrs['class'] = att.__class__.__name__


# define attribute_register

_tmp = [
    ('coordinates',),
    ('numbers',),
    ('obasis',),
    ('grid',),
    ('wfn',),
    ('cache', 'olp', 'o'),
    ('cache', 'kin', 'o'),
    ('cache', 'na', 'o'),
    ('cache', 'er', 'o'),
    ('extra',),
    ('pseudo_numbers',),
    ('cell',),
]

attribute_register = {}
for record in _tmp:
    chk_field = CHKField(*record)
    name = '.'.join(record[:2])
    attribute_register[name] = chk_field
