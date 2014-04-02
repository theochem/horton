#!/usr/bin/env python

import h5py as h5
from glob import glob

for fn_h5 in glob('*.h5'):
    with h5.File(fn_h5) as f:
        if 'dm_beta' in f['wfn']:
            dm_full = f['wfn/dm_alpha/array'][:] + f['wfn/dm_beta/array'][:]
            dm_spin = f['wfn/dm_alpha/array'][:] - f['wfn/dm_beta/array'][:]
            g = f.create_group('dm_spin')
            g['array'] = dm_spin
            g.attrs['class'] = 'DenseOneBody'
        else:
            dm_full = f['wfn/dm_alpha/array'][:]*2
        g = f.create_group('dm_full')
        g['array'] = dm_full
        g.attrs['class'] = 'DenseOneBody'
        del f['wfn']
        if 'props' in f:
            for key, value in f['props'].iteritems():
                f[key] = f['props'][key]
            del f['props']
