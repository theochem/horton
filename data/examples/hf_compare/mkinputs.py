#!/usr/bin/env python

import os

from horton import *  # pylint: disable=wildcard-import,unused-wildcard-import


for z in range(1, 18):
    for charge in range(-1, 3):
        nel = z-charge
        if nel <= 0:
            continue
        if nel%2 == 0:
            mults = [1, 3, 5, 7, 9]
        else:
            mults = [2, 4, 6, 8]
        for mult in mults:
            if (nel-(mult-1))/2 < 0:
                continue
            dirname = '%03i_%s_Q%+i_M%i' % (z, periodic[z].symbol.rjust(2, '_'), charge, mult)
            if not os.path.isdir(dirname):
                os.mkdir(dirname)
            with open('%s/gaussian.com' % dirname, 'w') as f:
                print('%nproc=4', file=f)
                print('%chk=gaussian.chk', file=f)
                print('#p HF/aug-cc-pVDZ scf(tight,xqc) IOP(6/7=3)', file=f)
                print(file=f)
                print(dirname, file=f)
                print(file=f)
                print(charge, mult, file=f)
                print(periodic[z].symbol, 0.0, 0.0, 0.0, file=f)
                print(file=f)
