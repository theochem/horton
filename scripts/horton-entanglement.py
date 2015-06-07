#!/usr/bin/env python
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


import sys, numpy as np


def plot(i12ind1, i12ind2, i12val, orbinit, orbfinal, thresh, s1index, s1value):
    try:
        import matplotlib.pyplot as plt
        import matplotlib.lines as mlines
        from matplotlib.ticker import NullFormatter
        import math
    except ImportError:
        if log.do_warning:
            log.warn('Skipping plots because matplotlib was not found.')
        return

    norb = orbfinal-orbinit
    orbitals = np.arange(orbinit, orbfinal)
    theta = 2 * np.pi * (orbitals-orbinit)/(norb)
    r = 22*np.ones(norb,int)-3.00*((orbitals-orbinit)%3)

    plt.figure(figsize=(10,5))
    ax = plt.subplot(121, polar=True)
    ax.grid(False)
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)

    plt.plot(theta, r, 'o', markersize=12, alpha=0.2)
    for i in range(len(orbitals)):
        plt.annotate(
            i+1+orbinit,
            xy = (theta[i], r[i]), xytext = (0, 0),
            textcoords = 'offset points', ha = 'center', va = 'bottom',
            fontsize=8, fontweight='bold',
            )

    ax.yaxis.set_data_interval(0,22.5)
    ax.xaxis.set_major_formatter(NullFormatter())
    ax.yaxis.set_major_formatter(NullFormatter())
    legend = []
    for ind in range(len(i12val)):
        if i12val[ind] >= thresh:
            if i12val[ind] >= 0.0001 and i12val[ind] < 0.001:
                plt.plot([theta[i12ind1[ind]-orbinit], theta[i12ind2[ind]-orbinit]], [r[i12ind1[ind]-orbinit], r[i12ind2[ind]-orbinit]], ':', lw=2, color='orange')
            if i12val[ind] >= 0.001 and i12val[ind] < 0.01:
                plt.plot([theta[i12ind1[ind]-orbinit], theta[i12ind2[ind]-orbinit]], [r[i12ind1[ind]-orbinit], r[i12ind2[ind]-orbinit]], '-.', lw=2, color='g')
            if i12val[ind] >= 0.01 and i12val[ind] < 0.1:
                plt.plot([theta[i12ind1[ind]-orbinit], theta[i12ind2[ind]-orbinit]], [r[i12ind1[ind]-orbinit], r[i12ind2[ind]-orbinit]], '--', lw=2, color='r')
            if i12val[ind] >= 0.1:
                plt.plot([theta[i12ind1[ind]-orbinit], theta[i12ind2[ind]-orbinit]], [r[i12ind1[ind]-orbinit], r[i12ind2[ind]-orbinit]], '-', lw=3, color='b')

    blue_line = mlines.Line2D([], [], color='blue', marker='', lw=3, ls='-', label='0.1')
    red_line = mlines.Line2D([], [], color='red', marker='', lw=2, ls='--', label='0.01')
    green_line = mlines.Line2D([], [], color='green', marker='', ls='-.', lw=2, label='0.001')
    orange_line = mlines.Line2D([], [], color='orange', marker='', ls=':', lw=2, label='0.0001')

    if thresh >= 0.0001 and thresh < 0.001:
        legend.append(blue_line)
        legend.append(red_line)
        legend.append(green_line)
        legend.append(orange_line)
    if thresh >= 0.001 and thresh < 0.01:
        legend.append(blue_line)
        legend.append(red_line)
        legend.append(green_line)
    if thresh >= 0.01 and thresh < 0.1:
        legend.append(blue_line)
        legend.append(red_line)
    if thresh >= 0.1:
        legend.append(blue_line)

    plt.legend(handles=legend, loc='center', fancybox=True, fontsize=10)
    plt.title('Mutual information')

    ax2 = plt.subplot(122)
    ax2.axis([orbinit, orbfinal, 0, 0.71])
    ax2.vlines(s1index, [0], s1value, color='r', linewidth=2, linestyle='-')
    plt.ylabel('single-orbital entropy')
    plt.xlabel('Orbital index')
    plt.plot(s1index, s1value, 'ro', markersize=8)

    plt.savefig('orbitalentanglement.png', dpi=300)
#   return plt

def read_i12_data(orbinit, orbfinal, thresh):
    index1 = np.array([])
    index2 = np.array([])
    value = np.array([])
    with open("i12.dat") as f:
        counter = 1
        for line in f:
            words = line.split()
            if len(words) != 3:
                raise IOError('Expecting 3 fields on each data line in i12.dat')
            if float(words[2]) >= thresh and int(words[0]) >= orbinit and int(words[1]) >= orbinit and int(words[0]) <= orbfinal and int(words[1]) <= orbfinal:
                index1 = np.append(index1, int(words[0])-1)
                index2 = np.append(index2, int(words[1])-1)
                value = np.append(value, float(words[2]))
    return index1, index2, value

def read_s1_data(orbinit, orbfinal):
    index = np.array([])
    value = np.array([])
    with open("s1.dat") as f:
        for line in f:
            words = line.split()
            if len(words) != 2:
                raise IOError('Expecting 2 fields on each data line in s1.dat')
            index = np.append(index, int(words[0]))
            value = np.append(value, float(words[1]))
    if orbfinal:
        newind = np.where(index>=(orbinit+1))
        index = index[newind]
        value = value[newind]
        newind2 = np.where(index<=orbfinal)
        index = index[newind2]
        value = value[newind2]
    return index, value

def main():

    args = sys.argv[1:]
    if len(args) == 0:
        print >> sys.stderr, 'Expecting at least one argument: threshold (float)'
        sys.exit(-1)
    thresh = float(args.pop(0))

    # Get plot range if provided, otherwise use full range
    try:
        orbinit = int(args.pop(0))-1
        orbfinal = int(args.pop(0))
    except:
        orbinit = 0
        orbfinal = None

    # Read s1.dat and store data
    s1index, s1value = read_s1_data(orbinit, orbfinal)

    if orbfinal is None:
        orbfinal = len(s1index)

    # Read i12.dat and store data
    i12index1, i12index2, i12value = read_i12_data(orbinit, orbfinal, thresh)

    # Plot i12 graph
    plt1 = plot(i12index1, i12index2, i12value, orbinit, orbfinal, thresh, s1index, s1value)
#   plt1.show()

if __name__ == '__main__':
    main()
