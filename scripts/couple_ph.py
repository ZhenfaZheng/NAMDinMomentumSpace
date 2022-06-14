#!/usr/bin/env python

import math, os
import numpy as np
import postnamd as pn
import matplotlib as mpl; mpl.use('agg')
import matplotlib.pyplot as plt
from glob import glob

def main():

    which_bas = [1]
    which_bas = range(300)

    qplabels = 'gkmg'
    qpath = np.array([
        [0.00000, 0.00000, 0.00000],
        [0.33333, 0.33333, 0.00000],
        [0.50000, 0.00000, 0.00000],
        [0.00000, 0.00000, 0.00000] ])

    inp = pn.read_inp('inp')

    prefix = inp['EPMPREF']; epmdir = inp['EPMDIR']
    filepm = os.path.join(epmdir, prefix + '_ephmat_p1.h5')
    qpts = pn.read_ephmath5(filepm, dset='/el_ph_band_info/q_list')
    phen = pn.read_ephmath5(filepm, dset='/el_ph_band_info/ph_disp_meV')

    coup = pn.read_couple(filcoup='EPECTXT', inp=inp)
    kkqmap = np.loadtxt('KKQMAP', dtype='int32')
    kkqmap -= 1

    nqs = qpts.shape[0]
    nmodes = coup.shape[0]
    coup_ph = np.zeros((nqs, nmodes))

    for ib in which_bas:
        coup_ph[kkqmap[ib], :] += coup[:, ib, :].T

    A = pn.read_ephmath5(filepm,
            dset='/el_ph_band_info/lattice_vec_angstrom')
    a1, a2, a3 = (A[0], A[1], A[2])
    b1, b2, b3 = pn.calc_rec_vec(a1, a2, a3)
    qpts_cart = pn.frac2cart(qpts, b1, b2, b3)
    qpath_cart = pn.frac2cart(qpath, b1, b2, b3)
    q_index, qp_index = pn.select_kpts_on_path(qpts, qpath, norm=0.001)
    q_loc, qp_loc = pn.loc_on_kpath(qpts_cart, q_index, qp_index, qpath_cart)

    plot_coup_ph(q_loc, phen, qp_loc, qplabels, coup_ph, q_index)


def plot_coup_ph(q_loc, phen, qp_loc, qplabels, coup_ph, index,
        figname='COUPLEPH.png'):

    X = q_loc
    E = phen[index, :]
    coup = coup_ph[index, :]
    nmodes = phen.shape[1]

    xmin = qp_loc[0] ; xmax = qp_loc[-1]
    ymin = E.min() ; ymax = E.max() ; dy = ymax - ymin
    ymin -= dy*0.05 ; ymax += dy*0.05

    figsize_x = 4.2
    figsize_y = 4.8 # in inches
    fig, ax = plt.subplots()
    fig.set_size_inches(figsize_x, figsize_y)
    mpl.rcParams['axes.unicode_minus'] = False

    cmap = 'plasma'
    cmin = coup.min() ; cmax = coup.max()
    norm = mpl.colors.Normalize(cmin, cmax)

    sort = np.argsort(X)
    for im in range(nmodes):
        ax.plot(X[sort], E[sort, im], '#1A5599', lw=0.7)

    nqpath = qp_loc.shape[0]
    for ipath in range(1, nqpath):
        x = qp_loc[ipath]
        ax.plot([x,x], [ymin,ymax], 'gray', lw=0.7, ls='--')

    for im in range(nmodes):
        sort = np.argsort(coup[:,im])
        sc = ax.scatter(X[sort], E[sort,im], s=10, lw=0,
                c=coup[sort,im], cmap=cmap, norm=norm)

    ticks = []
    for s in qplabels:
        s = u'\u0393' if (s=='g' or s=='G') else s.upper()
        ticks.append(s)

    ax.set_xticks(qp_loc)
    ax.set_xticklabels(ticks)

    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)

    ax.set_ylabel('Phonon energy (meV)')

    cbar = plt.colorbar(sc)
    cbar.set_label('Coupling (meV)')

    plt.tight_layout()
    plt.savefig(figname, dpi=400)
    print("\n%s has been saved.\n"%figname)


if __name__=='__main__':
    main()

