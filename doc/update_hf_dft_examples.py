#!/usr/bin/env python
'''Update the table with HF/DFT examples based on files in data/examples/hf_dft.'''

from glob import glob
from cStringIO import StringIO
import json, os
from common import write_rst_table, write_if_changed

cases = []
for fn_py in sorted(glob('../data/examples/hf_dft/*.py')):
    with open(fn_py) as f:
        s = ''
        for line in f:
            if line.startswith('#JSON'):
                s += line[5:]
    if len(s) == 0:
        raise RuntimeError('Could not find JSON line in HF/DFT example script.')
    meta = json.loads(s)
    case = [meta['difficulty'], os.path.basename(fn_py), meta['lot'],
            meta['scf'], meta['description']]
    cases.append(case)

cases.sort()
table = [['File name', 'LOT', 'SCF', 'Description']] + [case[1:] for case in cases]

f = StringIO()
write_rst_table(f, table)
write_if_changed('hf_dft_examples.rst.inc', f.getvalue())
