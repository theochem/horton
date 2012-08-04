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
"""This field defines all the parts of the System class than can be written
   to a checkpoint file, and how that should be done in each case.

   The philosophy of this module is to take away as much of the checkpointing
   machinery out of the rest of Horton. The only (somewhat unavoidable)
   exceptions to this philosophy are the ``load_checkpoint`` function in
   ``horton.io`` and the ``from_hdf5`` and ``to_hdf5`` methods in several parts
   of Horton.

   TODO: (long future) make the fields clever enough such that they can convert
   data from one linalgfactory into another.
"""


import h5py as h5, numpy as np



__all__ = ['register']


class CHKField(object):
    def __init__(self, att_name, key=None, att_class=None):
        """
           **Argument:**

           att_name
                The name of the attribute of this field

           **Optional argument:**

           key
                If the attribute is a dictionary, the key refers to an element
                in the dictionary. This must be a string

           att_class
                A class to create the attribute with when it is read from the
                checkpoint file. If this is an abstract class, the from_hdf5
                method of the abstract class must be clever enough to construct
                the right derived class. By default the name of the subclass is
                stored in the 'class' attribute of the corresponding group in
                the HDF5 file.

           The attribute of the System class may be None, a numpy array or
           an object that has to_hdf5 and from hdf5 methods method. In the last
           case, the att_classes argument must be given.
        """
        if not isinstance(att_name, basestring):
            raise TypeError('att_name must be a string')
        if key is not None and not isinstance(key, basestring):
            raise TypeError('key must be a string, when given')
        if att_class is not None and not isinstance(att_class, type):
            raise TypeError('att_class must be a class, when given.')
        self.att_name = att_name
        self.key = key
        self.att_class = att_class

    def read(self, chk):
        """Read and return an the attribute from the chk file

           **Argument:**

           chk
                A h5.File object corresponding to a Horton checkpoint file.

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
        if isinstance(item, h5.Dataset):
            # convert it to a numpy array
            return np.array(item)
        elif isinstance(item, h5.Group):
            if self.att_class is None:
                raise ValueError('The field matches a Group object in the checkpoint file but not att_class is given.')
            return att_class.from_hdf5(item)


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
            att = att.get(self.key)
            if att is None:
                return
        # B) Optionally make a datagroup
        grp = chk
        name = self.att_name
        if self.key is not None:
            grp = chk.require_group(self.att_name)
            name = self.key
        # C) Dump the data to HDF5
        if isinstance(att, np.ndarray):
            dset = grp.require_dataset(name, shape=att.shape, dtype=att.dtype)
            dset.write_direct(att)
        else:
            grp = grp.require_group(name)
            att.to_hdf5(grp)
            # needed to create object of the right type when reading from
            # checkpoint:
            grp.attrs['class'] = att.__class__.__name__


register = {
    'coordinates': CHKField('coordinates'),
    'numbers': CHKField('numbers'),
}
