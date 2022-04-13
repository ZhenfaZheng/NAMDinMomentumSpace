#!/usr/bin/env python

import math, os
import numpy as np
import postnamd as pn
import matplotlib as mpl; mpl.use('agg')
import matplotlib.pyplot as plt
from glob import glob

def main():

    Eref = 0.0
    inp = pn.read_inp('inp')

    coup = pn.read_couple(filcoup='EPECTXT', inp=inp)
    coup_av = np.average(np.abs(coup), axis=0)
    plot_couple(coup_av, figname='COUPLE.png')

    filshps = glob('SHPROP.*')
    shp = pn.readshp(filshps)
    en, kpts = pn.ek_selected(inp=inp)
    plot_tdprop(shp, Eref, lplot=2, ksen=en, figname='TDEN.png')

    prefix = inp['EPMPREF']; epmdir = inp['EPMDIR']
    filepm = os.path.join(epmdir, prefix + '_ephmat_p1.h5')
    A = pn.read_ephmath5(filepm, dset='/el_ph_band_info/lattice_vec_angstrom')
    a1, a2, a3 = (A[0], A[1], A[2])
    b1, b2, b3 = pn.calc_rec_vec(a1, a2, a3)
    kpts_cart = pn.frac2cart(kpts, b1, b2, b3)

    plot_kprop(kpts_cart, shp, B=[b1, b2, b3], axis='xy', figname='TDKPROPxy.png')
    plot_kprop(kpts_cart, shp, B=[b1, b2, b3], axis='yz', figname='TDKPROPyz.png')
    plot_kprop(kpts_cart, shp, B=[b1, b2, b3], axis='xz', figname='TDKPROPxz.png')


    kplabels = 'gkmg'
    kpath = np.array([
        [0.00000, 0.00000, 0.00000],
        [0.33333, 0.33333, 0.00000],
        [0.50000, 0.00000, 0.00000],
        [0.00000, 0.00000, 0.00000]
        ])
    kpath_cart = pn.frac2cart(kpath, b1, b2, b3)
    k_index, kp_index = pn.select_kpts_on_path(kpts, kpath, norm=0.01)
    k_loc, kp_loc = pn.loc_on_kpath(kpts_cart, k_index, kp_index, kpath_cart)

    plot_tdband(k_loc, en, kp_loc, kplabels, shp, k_index,
                Eref=Eref, figname='TDBAND.png')


def plot_tdband(k_loc, en, kp_loc, kplabels, shp, index,
        X_bg=None, E_bg=None, Eref=0.0, figname='TDBAND.png'):

    nbasis = index.shape[0]
    ntsteps = shp.shape[0]
    namdtime = shp[-1, 0]
    potim = namdtime / ntsteps
    pop = shp[:, index+2]

    X = np.tile(k_loc, ntsteps).reshape(ntsteps, nbasis)
    E = np.tile(en[index], ntsteps).reshape(ntsteps, nbasis) - Eref

    xmin = kp_loc[0] ; xmax = kp_loc[-1]
    ymin = E.min() ; ymax = E.max() ; dy = ymax - ymin
    ymin -= dy*0.05 ; ymax += dy*0.05

    figsize_x = 4.8
    figsize_y = 3.6 # in inches
    fig, ax = plt.subplots()
    fig.set_size_inches(figsize_x, figsize_y)
    mpl.rcParams['axes.unicode_minus'] = False

    cmap = 'rainbow'
    norm = mpl.colors.Normalize(0,namdtime)
    color_t = np.tile(np.arange(ntsteps), nbasis).reshape(nbasis,ntsteps).T * potim
    s_avg = np.average(pop[pop>0])
    dotsize = pop / s_avg * 5

    nkpath = kp_loc.shape[0]
    for ipath in range(1, nkpath):
        x = kp_loc[ipath]
        ax.plot([x,x], [ymin,ymax], 'k', lw=0.7, ls='--')

    sc = ax.scatter(X, E, s=dotsize, lw=0, c=color_t, cmap=cmap, norm=norm)
    cbar = plt.colorbar(sc, fraction=0.05)
    cbar.set_label('Time (fs)')

    ticks = []
    for s in kplabels:
        s = u'\u0393' if (s=='g' or s=='G') else s.upper()
        ticks.append(s)

    ax.set_xticks(kp_loc)
    ax.set_xticklabels(ticks)

    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.set_ylabel('Energy (eV)')

    plt.tight_layout()
    plt.savefig(figname, dpi=400)


def plot_couple(coup, figname='COUPLE.png'):
    '''
    This function plots average couplings.

    Parameters:
    coup: ndarray, average coupling data in forms of coup[nb, nb]
    figname: string, output figure file name.
    '''

    fig = plt.figure()
    figsize_x = 4.8
    figsize_y = 3.6 # in inches
    fig.set_size_inches(figsize_x, figsize_y)

    cmap = 'bwr'
    n = coup.shape[0]
    coup *= 1000.0 # change unit to meV
    Bmin = 0.5; Bmax = n + 0.5
    cmin = 0.0; cmax = np.max(coup)
    norm = mpl.colors.Normalize(cmin,cmax)
    plt.imshow(coup, cmap=cmap, origin='lower', norm=norm,
        extent=(Bmin,Bmax,Bmin,Bmax), interpolation='none')

    cbar = plt.colorbar()
    # cbar.ax.set_title('   meV')
    cbar.set_label('Coupling (meV)')
    plt.tight_layout()
    plt.savefig(figname, dpi=400)


def plot_tdprop(shp, Eref=0.0, lplot=1, ksen=None, figname='tdshp.png'):
    '''
    This function loads data from SHPROP.xxx files,
    and plot average evolution of energy & surface hopping proportion of
    electronic states.

    Parameters:
    shp    : ndarray, average data of SHPROP.xxx files, in forms of
             shp[ntsteps, nb+2].
    Eref   : float, energy reference. Make sure shp & ksen have same Eref!!!
    lplot  : integer, fig type to plot. 1: plot average proportions; 2: plot
             average energy evolution with proportions.
    ksen   : ndarray, KS energies in forms of ksen[nbands]. Here we suppose
             ksen do not change by time.
    figname: string, file name of output figure.
    '''

    figsize_x = 4.8
    figsize_y = 3.2 # in inches
    namdtime = shp[-1,0]

    fig, ax = plt.subplots()
    fig.set_size_inches(figsize_x, figsize_y)
    mpl.rcParams['axes.unicode_minus'] = False

    if lplot==1:
        ylabel = 'SHPROP'
        ax.plot(shp[:,0], shp[:,2:])
        ax.set_ylim(0.0,1.0)
        ax.set_ylabel(ylabel)
    else:
        cmap = 'hot_r'
        dotsize = 50
        ylabel = 'Energy (eV)'
        ntsteps = shp.shape[0]
        nbands = shp.shape[1] -2
        cmin = 0.0; cmax = np.max(shp[:,2:])
        cmax = math.ceil(cmax*10)/10
        norm = mpl.colors.Normalize(cmin,cmax)

        if (ksen.shape[0]!=nbands):
            print('\nNumber of ksen states doesn\'t match with SHPROP data!\n')
        E = np.tile(ksen-Eref, ntsteps).reshape(ntsteps, nbands)
        T = np.tile(shp[:,0], nbands).reshape(nbands,ntsteps).T
        sc = ax.scatter(T, E, s=dotsize, c=shp[:,2:], lw=0,
                        norm=norm, cmap=cmap)
        ax.plot(shp[:,1]-Eref, 'r', lw=1, label='Average Energy')
        plt.colorbar(sc)

        x1 = 0.05 * namdtime; x2 = 0.1 * namdtime
        for ib in range(nbands):
            y = ksen[ib] - Eref
            ax.plot([x1, x2], [y, y], color='r', lw=0.7)

        ax.set_ylabel(ylabel)

    ax.set_xlim(0,namdtime)
    ax.set_xlabel('Time (fs)')

    plt.tight_layout()
    plt.savefig(figname, dpi=400)


def plot_kprop(kpts, shp, B, axis='xy', figname='TDKPROP.png'):

    axdict = {'x':0, 'y':1, 'z':2}
    # axis must be set as 'xy', 'yz' or 'xz'!
    kax = [axdict[s] for s in axis]

    gamma = [0, 0, 0]
    Bpath = np.array([gamma, B[kax[0]], B[kax[0]]+B[kax[1]], B[kax[1]], gamma])
    xmin = Bpath[:,kax[0]].min(); xmax = Bpath[:,kax[0]].max()
    ymin = Bpath[:,kax[1]].min(); ymax = Bpath[:,kax[1]].max()

    namdtime = shp[-1,0]
    ntsteps = shp.shape[0]
    potim = namdtime / ntsteps
    nbasis = shp.shape[1] - 2

    scale = np.sqrt( 12.0 / ((ymax - ymin) * (xmax - xmin)) )
    figsize_x = (xmax - xmin) * scale * 1.3
    figsize_y = (ymax - ymin) * scale # in inches

    fig, ax = plt.subplots()
    fig.set_size_inches(figsize_x, figsize_y)
    mpl.rcParams['axes.unicode_minus'] = False

    cmap = 'rainbow'
    norm = mpl.colors.Normalize(0,namdtime)
    color_t = np.tile(np.arange(ntsteps), nbasis).reshape(nbasis,ntsteps).T * potim
    s_avg = np.average(shp[:,2:][shp[:,2:]>0])
    dotsize = shp[:,2:] / s_avg * 0.5

    X = np.tile(kpts[:,kax[0]], ntsteps).reshape(ntsteps,nbasis)
    Y = np.tile(kpts[:,kax[1]], ntsteps).reshape(ntsteps,nbasis)

    ax.plot(Bpath[:, kax[0]], Bpath[:, kax[1]], color='k', lw=1.0)

    sc = ax.scatter(X, Y, s=dotsize, lw=0, c=color_t, norm=norm, cmap=cmap)
    cbar = plt.colorbar(sc, fraction=0.05)
    cbar.set_label('Time (fs)')

    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.set_aspect('equal')

    plt.axis('off')
    plt.tight_layout()
    plt.savefig(figname, dpi=600)

if __name__=='__main__':
    main()
