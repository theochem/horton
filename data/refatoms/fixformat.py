#!/usr/bin/env python

import h5py as h5
from glob import glob

for fn_h5 in glob('*.h5'):
    with h5.File(fn_h5, "r+") as f:
        if 'dm_beta' in f['wfn']:
            dm_full = f['wfn/dm_alpha/array'][:] + f['wfn/dm_beta/array'][:]
            dm_spin = f['wfn/dm_alpha/array'][:] - f['wfn/dm_beta/array'][:]
            f['dm_spin'] = dm_spin
        else:
            dm_full = f['wfn/dm_alpha/array'][:]*2
        f['dm_full'] = dm_full
        del f['wfn']
        if 'props' in f:
            for key, value in f['props'].items():
                f[key] = f['props'][key]
            del f['props']
