#!/usr/bin/env python

import h5py as h5
from glob import glob

for fn_h5 in glob('*.h5'):
    with h5.File(fn_h5) as f:
        if f['wfn'].attrs['class'] == 'OpenShellWFN':
            f['wfn'].attrs['class'] = 'UnrestrictedWFN'
        elif f['wfn'].attrs['class'] == 'ClosedShellWFN':
            f['wfn'].attrs['class'] = 'RestrictedWFN'
