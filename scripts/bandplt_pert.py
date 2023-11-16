#!/usr/bin/env python

import numpy as np
import matplotlib as mpl; mpl.use('agg')
import matplotlib.pyplot as plt

# Parameter Setting
emin = -8
emax = 12
prefix = 'graphene'
Efermi = -4.5302
figname = 'Band_pert.png'
kpath = 'gmkg'
knum = [30, 30, 30, 1]


def main():

    nks = 0
    for ii, kk in enumerate(knum):
        knum[ii] = nks
        if kk==1:
            nks += 1
        else:
            nks += 1 + kk

    filbands = prefix + '.bands'
    bdata = np.loadtxt(filbands)
    kaxis = bdata[:nks, 0]
    eig = bdata[:,4].reshape(-1, nks).T
    eig -= Efermi

    fig, ax = plt.subplots()
    fig.set_size_inches(3.6,4.8)
    mpl.rcParams['axes.unicode_minus'] = False
    ax.plot(kaxis, eig, 'r', lw=0.7)
    for kk in knum:
        ax.axvline(kaxis[kk], emin, emax, c='k', ls='--', lw=0.7)
    ax.set_xlim(kaxis[0], kaxis[-1])
    ax.set_ylim(emin, emax)
    ax.set_xticks(kaxis[knum])
    ax.set_xticklabels(HSkptLabels(kpath))
    ax.set_ylabel('Energy (eV)')

    plt.tight_layout()
    plt.savefig(figname, dpi=250)


def HSkptLabels(kpath):
    '''
    Instruct k axis tick labels.
    '''
    kticklabels = []
    for s in kpath:
        if s=='g' or s=='G':
            s = u'\u0393'
        else:
            s = s.upper()
        kticklabels.append(s)
    return kticklabels


if __name__=='__main__':
    main()
