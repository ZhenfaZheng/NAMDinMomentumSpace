#!/usr/bin/env python

import math, os
import numpy as np
import postnamd as pn
import matplotlib as mpl; mpl.use('agg')
import matplotlib.pyplot as plt
from glob import glob

def main():

    #                       Basic parameters setting                      #
    #######################################################################

    Eref = 0.0
    which_plt = [1, 11, 12, 2, 31, 32, 33, 4, 5, 6, 7, 81]
    '''
    Select which figures to plot.
    1: COUPLE.png; 11:COUPLE_EL.png; 12:COUPLE_PH.png; 2: TDEN.png;
    31: TDKPROPxy.png; 32: TDKPROPyz.png; 33: TDKPROPxz.png;
    4: TDBAND.png; 5: TDPH.png; 6:TDPHEN.png; 7:TDPHNUM.png
    81: TDQPROPxy.png; 82: TDQPROPyz.png; 83: TDQPROPxz.png;
    9: DISTRIBUTION.png; 10: TPROP.png
    '''

    kplabels = 'gmkg'
    kpath = np.array([ # for TDBAND.png
        [0.00000, 0.00000, 0.00000],
        [0.50000, 0.00000, 0.00000],
        [0.33333, 0.33333, 0.00000],
        [0.00000, 0.00000, 0.00000]
        ])
    qpath = kpath # for TDPH.png
    qplabels = kplabels

    '''
    You can choose phonon modes to save their eigenvectors
    by setting ph_indices = [ [iq1, im1], [iq2, im2], ... ]
    '''
    # ph_indices = [[0, 0], [0,1], [0,2]]
    ph_indices = []
    '''
    If you set ph_indices = [], this script will pick primary phonon modes
    with larger couplings to save automatically.
    '''


    #                              Read data                              #
    #######################################################################

    print("\nReading data ...")
    inp = pn.read_inp('inp')
    prefix = inp['EPMPREF']; epmdir = inp['EPMDIR']
    filepm = os.path.join(epmdir, prefix + '_ephmat_p1.h5')
    A = pn.read_ephmath5(filepm, dset='/el_ph_band_info/lattice_vec_angstrom')

    en, kpts = pn.ek_selected(inp=inp) # en & kpts selected
    Enk = pn.get_Enk_tot(kpath, A, inp=inp) # total Enk
    qpts = pn.read_ephmath5(filepm, dset='/el_ph_band_info/q_list')
    phen = pn.read_ephmath5(filepm, dset='/el_ph_band_info/ph_disp_meV')

    k_index, k_loc, kp_loc = pn.pick_kpts_on_path(kpts, kpath, A, norm=0.001)
    q_index, q_loc, qp_loc = pn.pick_kpts_on_path(qpts, qpath, A, norm=0.001)

    if ((1 in which_plt) or (11 in which_plt)):

        if os.path.isfile('EPECTXT'):
            coup = pn.read_couple(filcoup='EPECTXT', inp=inp)
            coup = coup * 1000.0 # change unit to meV
            coup_av = np.average(np.abs(coup), axis=0)
        elif os.path.isfile('EPELTXT'):
            coup_av = np.loadtxt('EPELTXT') * 1000.0
        else:
            print("\nERROR: EPELTXT file is not found!")

    if ((12 in which_plt) or (len(ph_indices)==0)):
        if os.path.isfile('EPPHTXT'):
            coup_ph = np.loadtxt('EPPHTXT')
            coup_ph *= 1000.0 # change unit to meV
        else:
            print("\nERROR: EPPHTXT file is not found!")

    l_read_shp = False
    for ii in which_plt:
        if (ii in [2, 31, 32, 33, 4, 9, 10]):
            l_read_shp = True; break
    if (l_read_shp):
        filshps = glob('SHPROP.*')
        if filshps:
            shp = pn.readshp(filshps)
        else:
            shp = None
            print('\nERROR: SHPROP files are not found!')
        print('SHPROP files have been read!')

    l_read_php = False
    for ii in which_plt:
        if (ii in [5, 6, 7, 81, 82, 83]):
            l_read_php = True; break
    if (l_read_php):
        nmodes = phen.shape[1] ; nqs = phen.shape[0]
        filphps = glob('PHPROP.*')
        if filphps:
            filphps = ['PHPROP.%d'%(im+1) for im in range(nmodes)]
            php = pn.readphp(filphps)
            print('PHPROP files have been read!')
        elif os.path.isfile('PHPROP'):
            php = np.loadtxt('PHPROP')
            ntsteps = php.shape[0] / nmodes
            php = php.reshape(nmodes, ntsteps, nqs+2)
            print('PHPROP files have been read!')
        else:
            php = None
            print('\nERROR: PHPROP files are not found!')


    #                             Plot figures                            #
    #######################################################################

    if (1 in which_plt):
        plot_couple(coup_av, figname='COUPLE.png')
    if (11 in which_plt):
        plot_couple_el(coup_av, k_loc, en, kp_loc, kplabels, k_index, Enk,
                       Eref, figname='COUPLE_EL.png')
    if (12 in which_plt):
        plot_coup_ph(coup_ph, q_loc, phen, qp_loc, qplabels, q_index,
                     figname='COUPLE_PH.png')

    if (2 in which_plt):
        plot_tdprop(shp, Eref, lplot=2, ksen=en, figname='TDEN.png')

    # times = [0, 50, 100, 200, 500, 1000]
    namdtime = int( inp['NAMDTIME'] )
    times = list( range(0, namdtime+1, int(namdtime/5)) )

    if (31 in which_plt):
        plot_tdkprop(kpts, shp, times, axis='xy', figname='TDKPROPxy.png')
    if (32 in which_plt):
        plot_tdkprop(kpts, shp, times, axis='yz', figname='TDKPROPyz.png')
    if (33 in which_plt):
        plot_tdkprop(kpts, shp, times, axis='xz', figname='TDKPROPxz.png')

    if (4 in which_plt):
        plot_tdband_sns(k_loc, en, kp_loc, kplabels, shp, k_index, times,
                        Enk=Enk, Eref=Eref, figname='TDBAND.png')

    if (5 in which_plt):
        plot_tdph_sns(q_loc, phen, qp_loc, qplabels, php, q_index, times,
                      figname='TDPH.png')

    if (6 in which_plt):
        plot_tdphen(php, figname='TDPHEN.png')

    if (7 in which_plt):
        plot_tdphnum(php, figname='TDPHNUM.png')

    if (81 in which_plt):
        print('')
        # for im in range(nmodes):
        #     figname = 'TDQPROPxy_Mode%d.png'%(im+1)
        #     plot_tdqprop(qpts, php, times, im, axis='xy', figname=figname)
        plot_tdqprop(qpts, php, times, axis='xy', figname='TDQPROPxy_tot.png')

    if (82 in which_plt):
        print('')
        # for im in range(nmodes):
        #     figname = 'TDQPROPyz_Mode%d.png'%(im+1)
        #     plot_tdqprop(qpts, php, times, im, axis='yz', figname=figname)
        plot_tdqprop(qpts, php, times, axis='yz', figname='TDQPROPyz_tot.png')

    if (83 in which_plt):
        print('')
        # for im in range(nmodes):
        #     figname = 'TDQPROPxz_Mode%d.png'%(im+1)
        #     plot_tdqprop(qpts, php, times, im, axis='xz', figname=figname)
        plot_tdqprop(qpts, php, times, axis='xz', figname='TDQPROPxz_tot.png')

    if (9 in which_plt):
        plot_distrib(en, shp, times, Ef=Eref, figname='DISTRIBUTION.png')

    if (10 in which_plt):
        plot_Tprop(en, shp, Ef=Eref, figname='TPROP.png')

    if (len(ph_indices)==0):
        ph_indices = pick_prim_phmod(coup_ph, q_index)
        print("\nSave primary phonon modes with larger coupling.")
        save_phmod(inp, ph_indices, coup_ph)
    else:
        print("\nSave choosen phonon modes.")
        save_phmod(inp, ph_indices)

    print("\nDone!\n")


###############################################################################

def plot_couple_el(coup_in, k_loc, en, kp_loc, kplabels, index, Enk,
                   Eref=0.0, figname='COUPLE_EL.png'):

    X = k_loc
    E = en[index] - Eref
    coup = np.sum(coup_in, axis=1)[index]
    nbands = Enk.shape[1] - 1

    xmin = kp_loc[0] ; xmax = kp_loc[-1]
    ymin = E.min() ; ymax = E.max() ; dy = ymax - ymin
    ymin -= dy*0.05 ; ymax += dy*0.05

    # cmap = 'hot_r'
    # cmin = np.min(coup[coup>0.0]); cmax = np.max(coup)
    # norm = mpl.colors.LogNorm(cmin,cmax)
    cmap = 'rainbow'
    cmin = 0.0; cmax = np.max(coup)
    norm = mpl.colors.Normalize(cmin,cmax)

    figsize_x = 4.2
    figsize_y = 4.8 # in inches
    fig, ax = plt.subplots()
    fig.set_size_inches(figsize_x, figsize_y)
    mpl.rcParams['axes.unicode_minus'] = False

    sc = ax.scatter(X, E, s=30, lw=0,
            c=coup, cmap=cmap, norm=norm)

    for ib in range(nbands):
        ax.plot(Enk[:,0], Enk[:,ib+1]-Eref, '#1A5599', lw=0.7)

    nkpath = kp_loc.shape[0]
    for ipath in range(1, nkpath):
        x = kp_loc[ipath]
        ax.plot([x,x], [ymin,ymax], 'k', lw=0.7, ls='--')

    ticks = []
    for s in kplabels:
        s = u'\u0393' if (s=='g' or s=='G') else s.upper()
        ticks.append(s)

    ax.set_xticks(kp_loc)
    ax.set_xticklabels(ticks)

    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.set_ylabel('Energy (eV)')

    cbar = plt.colorbar(sc, aspect=30)
    cbar.set_label('Coupling (meV)')
    plt.tight_layout()
    plt.savefig(figname, dpi=400)
    plt.close(fig)
    print("\n%s has been saved."%figname)


def plot_couple(coup, figname='COUPLE.png'):
    '''
    This function plots average couplings.

    Parameters:
    coup   : ndarray, average coupling data in forms of coup[nb, nb]
    figname: string, output figure file name.
    '''

    fig = plt.figure()
    figsize_x = 4.8
    figsize_y = 3.6 # in inches
    fig.set_size_inches(figsize_x, figsize_y)

    cmap = 'bwr'
    n = coup.shape[0]
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
    plt.close(fig)
    print("\n%s has been saved."%figname)


def plot_coup_ph(coup_ph, q_loc, phen, qp_loc, qplabels, index,
        figname='COUPLE_PH.png'):

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

    cmap = 'rainbow'
    cmin = coup.min() ; cmax = coup.max()
    norm = mpl.colors.Normalize(cmin, cmax)
    # cmap = 'plasma'
    # cmin = np.min(coup[coup>0.0]); cmax = np.max(coup)
    # norm = mpl.colors.LogNorm(cmin,cmax)

    sort = np.argsort(X)
    for im in range(nmodes):
        ax.plot(X[sort], E[sort, im], '#1A5599', lw=0.7)

    nqpath = qp_loc.shape[0]
    for ipath in range(1, nqpath):
        x = qp_loc[ipath]
        ax.plot([x,x], [ymin,ymax], 'gray', lw=0.7, ls='--')

    for im in range(nmodes):
        sort = np.argsort(coup[:,im])
        sc = ax.scatter(X[sort], E[sort,im], s=30, lw=0,
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

    cbar = plt.colorbar(sc, aspect=30)
    cbar.set_label('Coupling (meV)')

    plt.tight_layout()
    plt.savefig(figname, dpi=400)
    plt.close(fig)
    print("\n%s has been saved."%figname)


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
        dotsize = 30
        ylabel = 'Energy (eV)'
        ntsteps = shp.shape[0]
        nbands = shp.shape[1] -2

        # cmin = 0.0; cmax = np.max(shp[:,2:])
        # cmax = math.ceil(cmax*10)/10
        # norm = mpl.colors.Normalize(cmin,cmax)
        pop = shp[:,2:]
        cmin = 1.0e-4; cmax=1.0
        # cmin = np.min(pop[pop>0.0]); cmax = np.max(pop)
        norm = mpl.colors.LogNorm(cmin,cmax)

        if (ksen.shape[0]!=nbands):
            print('\nNumber of ksen states doesn\'t match with SHPROP data!\n')
        sort = np.argsort(np.sum(pop, axis=0))
        E = np.tile(ksen-Eref, ntsteps).reshape(ntsteps, nbands)
        T = np.tile(shp[:,0], nbands).reshape(nbands,ntsteps).T
        sc = ax.scatter(T, E[:,sort], s=dotsize, c=pop[:, sort], lw=0,
                        norm=norm, cmap=cmap)
        ax.plot(shp[:,0], shp[:,1]-Eref, 'b', lw=1, label='Average Energy')
        plt.colorbar(sc)

        # x1 = 0.05 * namdtime; x2 = 0.1 * namdtime
        # for ib in range(nbands):
        #     y = ksen[ib] - Eref
        #     ax.plot([x1, x2], [y, y], color='r', lw=0.7)

        ax.set_ylabel(ylabel)

    ax.set_xlim(0,namdtime)
    ax.set_xlabel('Time (fs)')

    plt.tight_layout()
    plt.savefig(figname, dpi=400)
    plt.close(fig)
    print("\n%s has been saved."%figname)


def plot_tdkprop(kpts, shp, times, axis='xy', figname='TDKPROP.png'):

    axdict = {'x':0, 'y':1, 'z':2}
    # axis must be set as 'xy', 'yz' or 'xz'!
    kax = [axdict[s] for s in axis]

    tindex = times2index(times, shp)
    times = shp[tindex,0]
    nts = len(times)
    pop = shp[tindex,2:]
    cmin = np.min(pop[pop>0.0]); cmax = np.max(pop)
    norm = mpl.colors.LogNorm(cmin,cmax)
    cmap = 'hot_r'

    if (nts < 4):
        ncol = nts; nrow = 1
    else:
        ncol = math.ceil(np.sqrt(nts))
        nrow = math.ceil( nts / ncol )

    figsize_x = 2.7 * ncol + 1.0
    figsize_y = 2.8 * nrow # in inches
    fig, axes = plt.subplots(nrow, ncol)
    fig.set_size_inches(figsize_x, figsize_y)
    mpl.rcParams['axes.unicode_minus'] = False

    X = kpts[:,kax[0]]
    Y = kpts[:,kax[1]]

    if (nts < ncol * nrow):
        for it in range(nts, ncol*nrow):
            ax = axes[math.floor(it/ncol), it%ncol]
            ax.axis('off')

    for it in range(nts):

        if (nrow==1):
            ax = axes[it] if ncol>1 else axes
        else:
            ax = axes[math.floor(it/ncol), it%ncol]

        ax.set_title('%.0f fs'%times[it])

        sort = np.argsort(pop[it,:])
        sc = ax.scatter(X[sort], Y[sort], s=10, lw=0, c=pop[it,sort],
                        norm=norm, cmap=cmap)

        ax.set_xlim(-0.02, 1.01)
        ax.set_ylim(-0.02, 1.01)
        ax.set_aspect('equal')
        # ax.set_xticks([])
        # ax.set_yticks([])
        ax.set_xlabel('k$_%s$'%axis[0])
        ax.set_ylabel('k$_%s$'%axis[1])

    ll = 0.6/figsize_x ; rr = 1.0 - 0.7 / figsize_x
    bb = 0.5/figsize_y ; tt = 1.0 - 0.4 / figsize_y
    if (nts < ncol * nrow): rr = 1.0 - 0.2 / figsize_x
    fig.subplots_adjust(bottom=bb, top=tt, left=ll, right=rr,
                        wspace=0.4, hspace=0.45)
    w_cb = 0.12 / figsize_x
    l, b, w, h = ax.get_position().bounds
    cb_ax = fig.add_axes([l+w+0.4*w_cb, b, w_cb, h])
    cbar = plt.colorbar(sc, cax=cb_ax)

    plt.savefig(figname, dpi=400)
    plt.close(fig)
    print("\n%s has been saved."%figname)


def plot_tdband_sns(k_loc, en, kp_loc, kplabels, shp, index, times,
        b_index=None, Enk=None, Eref=0.0, figname='TDBAND.png'):

    tindex = times2index(times, shp)
    times = shp[tindex,0]
    pop = shp[:,index+2][tindex,:]
    nts = len(times)

    X = k_loc
    E = en[index] - Eref

    xmin = kp_loc[0] ; xmax = kp_loc[-1]
    ymin = E.min() ; ymax = E.max() ; dy = ymax - ymin
    ymin -= dy*0.05 ; ymax += dy*0.05

    if (nts < 4):
        ncol = nts; nrow = 1
    else:
        ncol = math.ceil(np.sqrt(nts))
        nrow = math.ceil( nts / ncol )

    figsize_x = 2.6 * ncol + 1.0
    figsize_y = 3.5 * nrow # in inches
    fig, axes = plt.subplots(nrow, ncol)
    fig.set_size_inches(figsize_x, figsize_y)
    mpl.rcParams['axes.unicode_minus'] = False

    cmap = 'hot_r'
    # cmin = np.min(pop[pop>0.0]); cmax = np.max(pop)
    cmin = 1.0e-4; cmax=1.0
    norm = mpl.colors.LogNorm(cmin,cmax)

    if (nts < ncol * nrow):
        for it in range(nts, ncol*nrow):
            ax = axes[math.floor(it/ncol), it%ncol]
            ax.axis('off')

    for it in range(nts):

        if (nrow==1):
            ax = axes[it] if ncol>1 else axes
        else:
            ax = axes[math.floor(it/ncol), it%ncol]

        ax.set_title('%.0f fs'%times[it])

        # plot background E-k lines
        if b_index is None:
            for ib in range(Enk.shape[1]-1):
                ax.plot(Enk[:,0], Enk[:,ib+1]-Eref, '#1A5599', lw=0.7)
        else:
            if (it==0):
                sort_index = np.argsort(k_loc)
                X_bg = k_loc[sort_index]
                E_bg = en[index][sort_index] - Eref
                b_ind_bg = b_index[index][sort_index]

            for ib in np.unique(b_ind_bg):
                ax.plot(X_bg[b_ind_bg==ib], E_bg[b_ind_bg==ib], '#1A5599', lw=0.7)

        sc = ax.scatter(X, E, s=50, lw=0,
                c=pop[it,:], cmap=cmap, norm=norm)

        nkpath = kp_loc.shape[0]
        for ipath in range(1, nkpath):
            x = kp_loc[ipath]
            ax.plot([x,x], [ymin,ymax], 'k', lw=0.7, ls='--')

        ticks = []
        for s in kplabels:
            s = u'\u0393' if (s=='g' or s=='G') else s.upper()
            ticks.append(s)

        ax.set_xticks(kp_loc)
        ax.set_xticklabels(ticks)

        ax.set_xlim(xmin, xmax)
        ax.set_ylim(ymin, ymax)
        if (it%ncol==0):
            ax.set_ylabel('Energy (eV)')

    ll = 0.8/figsize_x ; rr = 1.0 - 0.8 / figsize_x
    bb = 0.4/figsize_y ; tt = 1.0 - 0.4 / figsize_y
    if (nts < ncol * nrow): rr = 1.0 - 0.2 / figsize_x
    fig.subplots_adjust(bottom=bb, top=tt, left=ll, right=rr,
                        wspace=0.45, hspace=0.4)
    w_cb = 0.12 / figsize_x
    l, b, w, h = ax.get_position().bounds
    cb_ax = fig.add_axes([l+w+0.4*w_cb, b, w_cb, h])
    cbar = plt.colorbar(sc, cax=cb_ax)

    plt.savefig(figname, dpi=400)
    plt.close(fig)
    print("\n%s has been saved."%figname)


def plot_tdph_sns(q_loc, phen, qp_loc, qplabels, php, index, times,
                  X_bg=None, E_bg=None, figname='TDPH.png'):

    nmodes = php.shape[0]
    tindex = times2index(times, php[0,:,:])
    times = php[0,tindex,0]
    php = np.cumsum(php, axis=1)

    pop = php[:,:,index+2][:,tindex,:]
    nts = len(times)

    # calculate ph number variation between it-1 and it.
    for ii in range(nts-1, 1, -1):
        pop[:,ii,:] = pop[:,ii,:] - pop[:,ii-1,:]

    X = q_loc
    E = phen[index, :]

    xmin = qp_loc[0] ; xmax = qp_loc[-1]
    ymin = E.min() ; ymax = E.max() ; dy = ymax - ymin
    ymin -= dy*0.05 ; ymax += dy*0.05

    if (nts < 4):
        ncol = nts; nrow = 1
    else:
        ncol = math.ceil(np.sqrt(nts))
        nrow = math.ceil( nts / ncol )

    figsize_x = 2.6 * ncol + 1.0
    figsize_y = 3.5 * nrow # in inches
    fig, axes = plt.subplots(nrow, ncol)
    fig.set_size_inches(figsize_x, figsize_y)
    mpl.rcParams['axes.unicode_minus'] = False

    # cmin = 1.0e-4; cmax=1.0
    cmin = pop.min() ; cmax = pop.max()
    norm = mpl.colors.Normalize(cmin, cmax)
    cmap = genrcmap(cmin, cmax)
    # cmap = 'hot_r'
    # cmax = np.max(pop) * 1.2 ; cmin = cmax / 1000
    # norm = mpl.colors.LogNorm(cmin,cmax)
    size = np.sqrt(np.abs(pop))
    s_avg = np.average(size[size>0])
    size = size / size.max()  * 50

    if (nts < ncol * nrow):
        for it in range(nts, ncol*nrow):
            ax = axes[math.floor(it/ncol), it%ncol]
            ax.axis('off')

    for it in range(nts):

        if (nrow==1):
            ax = axes[it] if ncol>1 else axes
        else:
            ax = axes[math.floor(it/ncol), it%ncol]

        if (it==0):
            ax.set_title('%.0f fs'%times[it])
        else:
            ax.set_title('%.0f - %.0f fs'%(times[it-1], times[it]))

        sort_index = np.argsort(q_loc)
        for im in range(nmodes):
            # ax.plot(q_loc[sort_index], E[sort_index, im], '#1A5599', lw=0.7)
            ax.plot(q_loc[sort_index], E[sort_index, im], 'gray', lw=0.7)
        nqpath = qp_loc.shape[0]
        for ipath in range(1, nqpath):
            x = qp_loc[ipath]
            ax.plot([x,x], [ymin,ymax], 'gray', lw=0.7, ls='--')

        for im in range(nmodes):
            sort = np.argsort(pop[im,it,:])
            sc = ax.scatter(X[sort], E[sort,im], s=size[im,it,sort], lw=0,
            # sc = ax.scatter(X[sort], E[sort,im], s=50, lw=0,
                            c=pop[im,it,sort], cmap=cmap, norm=norm)

        ticks = []
        for s in qplabels:
            s = u'\u0393' if (s=='g' or s=='G') else s.upper()
            ticks.append(s)

        ax.set_xticks(qp_loc)
        ax.set_xticklabels(ticks)

        ax.set_xlim(xmin, xmax)
        ax.set_ylim(ymin, ymax)
        if(it%ncol==0): ax.set_ylabel('Phonon energy (meV)')

    ll = 0.8/figsize_x ; rr = 1.0 - 1.0 / figsize_x
    bb = 0.4/figsize_y ; tt = 1.0 - 0.4 / figsize_y
    if (nts < ncol * nrow): rr = 1.0 - 0.2 / figsize_x
    fig.subplots_adjust(bottom=bb, top=tt, left=ll, right=rr,
                        wspace=0.45, hspace=0.4)
    w_cb = 0.12 / figsize_x
    l, b, w, h = ax.get_position().bounds
    cb_ax = fig.add_axes([l+w+0.4*w_cb, b, w_cb, h])
    cbar = plt.colorbar(sc, cax=cb_ax)

    plt.savefig(figname, dpi=400)
    plt.close(fig)
    print("\n%s has been saved."%figname)


def genrcmap(cmin=-1.0, cmax=1.0):

    num_pos = 100
    num_neg = int( np.abs(cmin) / cmax * num_pos )
    index = np.hstack((np.linspace(0.1, 0.5, num_neg),
                       np.linspace(0.5, 0.9, num_pos)))
    if (cmin>=0): index = np.linspace(0.5, 0.9, num_pos)

    cmap = plt.cm.jet
    clist = [cmap(i) for i in index]

    cmap_new = mpl.colors.LinearSegmentedColormap.from_list(
                   name='mycmap', colors=clist)

    return cmap_new


def plot_tdphen(php, figname='TDPHEN.png'):

    phen = php[:, :, 1].T
    nmodes = phen.shape[1]
    X = php[0, :, 0]
    namdtime = X[-1]

    figsize_x = 4.8
    figsize_y = 3.2 # in inches

    fig, ax = plt.subplots()
    fig.set_size_inches(figsize_x, figsize_y)
    mpl.rcParams['axes.unicode_minus'] = False

    cmap = mpl.colormaps['rainbow']
    cls = [cmap(i) for i in np.linspace(0, 1, nmodes)]
    # cls = ['C%d'%im for im in range(nmodes)] # 'CN' colors
    lbs = ['mode %d'%(im+1) for im in range(nmodes)]

    phen = np.cumsum(np.cumsum(phen, axis=0), axis=1)
    ax.fill_between(X, phen[:, 0],
        color=cls[0], label=lbs[0], lw=0.3, alpha=0.75)
    for im in range(1, nmodes):
        ax.fill_between(X, phen[:, im], phen[:, im-1],
            color=cls[im], label=lbs[im], lw=0.3, alpha=0.75)

    if (nmodes<=100):
        ncol = int(np.sqrt(nmodes) / 2) + 1
        fsize = 10 - np.sqrt(nmodes) / 2
        ax.legend(loc=1, ncol=ncol, fontsize=fsize)

    ax.set_xlim(0,namdtime)
    ax.set_xlabel('Time (fs)')
    ax.set_ylabel('Energy (eV)')

    plt.tight_layout()
    plt.savefig(figname, dpi=400)
    plt.close(fig)
    print("\n%s has been saved."%figname)



def plot_tdphnum(php, figname='TDPHNUM.png'):

    phnum = np.sum(php[:,:,2:], axis=2).T
    nmodes = phnum.shape[1]
    X = php[0, :, 0]
    namdtime = X[-1]

    figsize_x = 4.8
    figsize_y = 3.2 # in inches

    fig, ax = plt.subplots()
    fig.set_size_inches(figsize_x, figsize_y)
    mpl.rcParams['axes.unicode_minus'] = False

    cmap = mpl.colormaps['rainbow']
    cls = [cmap(i) for i in np.linspace(0, 1, nmodes)]
    lbs = ['mode %d'%(im+1) for im in range(nmodes)]

    for im in range(nmodes):
        ax.plot(X, phnum[:,im], label=lbs[im], color=cls[im], lw=1.0)

    if (nmodes<=100):
        ncol = int(np.sqrt(nmodes) / 2) + 1
        fsize = 10 - np.sqrt(nmodes) / 2
        ax.legend(loc=1, ncol=ncol, fontsize=fsize)

    ax.set_xlim(0,namdtime)
    ax.set_xlabel('Time (fs)')
    ax.set_ylabel('Phonon Number')

    plt.tight_layout()
    plt.savefig(figname, dpi=400)
    plt.close(fig)
    print("\n%s has been saved."%figname)


def plot_tdqprop(qpts, php, times, im=None, axis='xy', figname='TDQPROP.png'):

    axdict = {'x':0, 'y':1, 'z':2}
    # axis must be set as 'xy', 'yz' or 'xz'!
    qax = [axdict[s] for s in axis]

    tindex = times2index(times, php[0,:,:])
    times = php[0,tindex,0]
    nts = len(times)

    if im is None:
        php = np.sum(np.cumsum(php, axis=1), axis=0)
        pop = php[tindex,2:]
    else:
        php = np.cumsum(php, axis=1)
        pop = php[im, tindex,2:]

    # calculate ph number variation between it-1 and it.
    for ii in range(nts-1, 1, -1):
        pop[ii,:] = pop[ii,:] - pop[ii-1,:]

    # cmin = np.min(pop[pop>0.0]); cmax = np.max(pop)
    # cmin = cmax * 1.0e-4
    # norm = mpl.colors.LogNorm(cmin,cmax)
    # cmap = 'hot_r'
    cmin = pop.min() ; cmax = pop.max()
    norm = mpl.colors.Normalize(cmin, cmax)
    cmap = genrcmap(cmin, cmax)
    size = np.sqrt(np.abs(pop))
    s_avg = np.average(size[size>0])
    size = size / size.max()  * 50

    if (nts < 4):
        ncol = nts; nrow = 1
    else:
        ncol = math.ceil(np.sqrt(nts))
        nrow = math.ceil( nts / ncol )

    figsize_x = 2.7 * ncol + 1.0
    figsize_y = 2.8 * nrow # in inches
    fig, axes = plt.subplots(nrow, ncol)
    fig.set_size_inches(figsize_x, figsize_y)
    mpl.rcParams['axes.unicode_minus'] = False

    X = qpts[:,qax[0]]
    Y = qpts[:,qax[1]]

    if (nts < ncol * nrow):
        for it in range(nts, ncol*nrow):
            ax = axes[math.floor(it/ncol), it%ncol]
            ax.axis('off')

    for it in range(nts):

        if (nrow==1):
            ax = axes[it] if ncol>1 else axes
        else:
            ax = axes[math.floor(it/ncol), it%ncol]

        if (it==0):
            ax.set_title('%.0f fs'%times[it])
        else:
            ax.set_title('%.0f - %.0f fs'%(times[it-1], times[it]))

        # sort = np.argsort(pop[it,:])
        # sc = ax.scatter(X[sort], Y[sort], s=10, lw=0, c=pop[it,sort],
        #                 norm=norm, cmap=cmap)
        sort = np.argsort(size[it,:])
        sc = ax.scatter(X[sort], Y[sort],
                 s=size[it,sort], lw=0, c=pop[it,sort],
                 norm=norm, cmap=cmap)

        ax.set_xlim(-0.02, 1.01)
        ax.set_ylim(-0.02, 1.01)
        ax.set_aspect('equal')
        # ax.set_xticks([])
        # ax.set_yticks([])
        ax.set_xlabel('q$_%s$'%axis[0])
        ax.set_ylabel('q$_%s$'%axis[1])

    ll = 0.6/figsize_x ; rr = 1.0 - 0.7 / figsize_x
    bb = 0.5/figsize_y ; tt = 1.0 - 0.4 / figsize_y
    if (nts < ncol * nrow): rr = 1.0 - 0.2 / figsize_x
    fig.subplots_adjust(bottom=bb, top=tt, left=ll, right=rr,
                        wspace=0.4, hspace=0.45)
    w_cb = 0.12 / figsize_x
    l, b, w, h = ax.get_position().bounds
    cb_ax = fig.add_axes([l+w+0.4*w_cb, b, w_cb, h])
    cbar = plt.colorbar(sc, cax=cb_ax)

    plt.savefig(figname, dpi=400)
    plt.close(fig)
    print("%s has been saved."%figname)


def plot_distrib(en, shp, times, Ef=0.0, figname='DISTRIBUTION.png'):

    sort = np.argsort(en)
    en = en[sort] - Ef
    enfit = np.arange(en.min(), en.max(), 0.01)

    tindex = times2index(times, shp)
    times = shp[tindex,0]
    pop = shp[tindex, 2:][:,sort]

    nsns = len(times)

    figsize_x = 4.8
    figsize_y = 1.2 * nsns # in inches

    fig, axes = plt.subplots(nsns, 1)
    fig.set_size_inches(figsize_x, figsize_y)
    fig.subplots_adjust(hspace=0.1)
    mpl.rcParams['axes.unicode_minus'] = False

    for  ii, time in enumerate(times):

        from scipy.optimize import curve_fit
        popt, pcov = curve_fit(pn.func_fd, en, pop[ii], p0=[1000])
        fdfit = pn.func_fd(enfit, *popt)

        ax = axes if (nsns==1) else axes[ii]

        ax.plot([0,0], [-0.1,1.1], 'k', ls='--', lw=0.5)
        ax.plot(enfit, fdfit, lw=1.0, color='b', alpha=0.8)
        ax.scatter(en, pop[ii], s=8.0, lw=0.0, color='r', alpha=0.5)

        s = "%d fs, %d K"%(time, popt[0])
        ax.text(0.97, 0.90, s, ha='right', va='top',
                bbox=dict(boxstyle='round', ec='k', alpha=0.1),
                transform=ax.transAxes
                )

        ax.set_ylim(-0.1, 1.1)
        ax.set_ylabel('Distribution')

        if ii<nsns-1:
            ax.set_xticks([])
        else:
            ax.set_xlabel('$E-E_f$ (eV)')

    plt.tight_layout()
    plt.savefig(figname, dpi=400)
    plt.close(fig)
    print("\n%s has been saved."%figname)


def plot_Tprop(en, shp, Ef=0.0, figname='TPROP.png'):

    en -= Ef
    ntsteps = shp.shape[0]

    from scipy.optimize import curve_fit
    itime = np.arange(0, ntsteps, 20)
    time = shp[itime, 0]
    pop = shp[:, 2:]

    temp = []
    for it in itime:
        popt, pcov = curve_fit(pn.func_fd, en, pop[it], p0=[1000])
        temp.append(popt[0])
    temp = np.array(temp)

    figsize_x = 4.8
    figsize_y = 3.2 # in inches

    fig, ax = plt.subplots()
    fig.set_size_inches(figsize_x, figsize_y)
    mpl.rcParams['axes.unicode_minus'] = False

    ax.plot(time, temp, lw=1.0)

    ax.set_xlabel('Time (fs)')
    ax.set_ylabel('Temperature (K)')

    plt.tight_layout()
    plt.savefig(figname, dpi=400)
    plt.close(fig)
    print("\n%s has been saved."%figname)


def times2index(times, shp):

    tindex = []
    for time in times:
        if (time==0 and shp[0,0]>0): time = shp[0,0]
        index = np.argwhere(shp[:,0]==time)
        if (len(index) > 0):
            tindex.append(index[0][0])
        else:
            print("\nWARNING: No data in SHPROP/PHPROP files " +\
                  "for time=%.0f fs!"%time)

    return np.array(tindex, dtype='int')


def pick_prim_phmod(coup_ph, q_index):

    # Choose primary phonon modes
    coup_max = coup_ph[q_index,:].max()
    coup_norm = coup_max / 10.0
    nmodes = coup_ph.shape[1]
    nqs = q_index.shape[0]
    nqparts = 1 if (nqs<=4) else 4
    nqs_p = int(nqs/nqparts)
    ph_indices = []
    for im in range(nmodes):
        coup = coup_ph[q_index, im]
        for ipart in range(nqparts):
            ist = ipart * nqs_p
            ied = (ipart+1)*nqs_p if (ipart<nqparts-1) else nqs
            coup_p = coup[ist:ied]
            ii = np.argpartition(-coup_p, 1)[0] + ist
            if (coup[ii] > coup_norm):
                iq = q_index[ii]
                ph_indices.append([iq,im])

    return ph_indices


def save_phmod(inp, ph_indices, coup_ph=None):

    ph_indices = np.array(ph_indices)
    nph = ph_indices.shape[0]
    if (nph == 0):
        print("\nNo phonon mode to save!!!")
        return

    prefix = inp['EPMPREF']; epmdir = inp['EPMDIR']
    filepm = os.path.join(epmdir, prefix + '_ephmat_p1.h5')
    qpts = pn.read_ephmath5(filepm, dset='/el_ph_band_info/q_list')
    phen = pn.read_ephmath5(filepm, dset='/el_ph_band_info/ph_disp_meV')
    latt_vec = pn.read_ephmath5(filepm, dset='/el_ph_band_info/lattice_vec_angstrom')
    phmod_ev_tot = pn.read_ephmath5(filepm, dset='/el_ph_band_info/phmod_ev_r')
    at_pos = pn.read_ephmath5(filepm, dset='/el_ph_band_info/atom_pos')
    at_mass = pn.read_ephmath5(filepm, dset='/el_ph_band_info/mass_a.u.')

    # Atom number maybe wrong, please check output .xsf files!!!
    at_num = [int(mass/2) for mass in at_mass]

    for iph in range(nph):

        iq = ph_indices[iph,0]
        im = ph_indices[iph,1]
        qpt = qpts[iq]
        phmod_ev = phmod_ev_tot[:,:,im,iq].T

        if (coup_ph is None):
            header  = "Phonon Mode\n"
            header += "# qpt: ( %10.6f%10.6f%10.6f ) iq: %d im: %d\n"%(*qpt, iq+1, im+1)
            header += "# freq: %.3f meV"%phen[iq,im]
        else:
            header  = "Phonon Mode with Larger Coupling\n"
            header += "# qpt: ( %10.6f%10.6f%10.6f ) iq: %d im: %d\n"%(*qpt, iq+1, im+1)
            header += "# freq: %.3f meV coupling: %.3G meV"%(phen[iq,im], coup_ph[iq,im])

        filname = 'PHMOD_%d_%d.xsf'%(iq+1, im+1)
        gen_xsf(latt_vec, at_pos, at_num, phmod_ev, header, filname)


def gen_xsf(latt_vec, at_pos, at_num, phmod_ev, header='', filname='PHMOD.xsf'):

    nat = at_pos.shape[0]
    at_pos_cart = np.matmul(at_pos, latt_vec)
    pos_phev = np.hstack((at_pos_cart, phmod_ev))

    mode = 'w'
    with open(filname, mode) as ff:

        if (header != ''):
            ff.write("# %s\n"%header)

        ff.write("CRYSTAL\n")
        ff.write("PRIMVEC\n")

        for vec in latt_vec:
            ff.write("  {:16.10f}{:16.10f}{:16.10f}\n".format(*vec))

        ff.write("PRIMCOORD\n")
        ff.write(" {:d} {:d}\n".format(nat, 1))

        coord_form = "{:d} {:12.6f}{:12.6f}{:12.6f}{:12.6f}{:12.6f}{:12.6f}\n"
        for iat in range(nat):
            ff.write(coord_form.format(at_num[iat], *pos_phev[iat]))

    print("%s has been saved."%filname)


if __name__=='__main__':
    main()
