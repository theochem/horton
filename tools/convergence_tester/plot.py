#!/usr/bin/env python

import matplotlib.pyplot as pt, numpy as np


with open('scf_results.txt') as f:
    pt.clf()
    tally  = dict()

    for line in f:
        words = line.split()
        mixing = float(words[0])

        temp = words[1]
        if temp == "inf":
            temp = 10
        niter = int(temp)

        if mixing not in tally:
            tally[mixing]=[]
        tally[mixing].append(niter)

    for k in sorted(tally.iterkeys()):
        tally[k].sort()
        pt.plot(tally[k], np.linspace(0.0, 1.0, len(tally[k])), label=str(-np.log10(k)), ls='steps-post')

    pt.legend(loc=0)
    pt.savefig('scf_results.png')
