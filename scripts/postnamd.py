#!/usr/bin/env python3

import h5py
import numpy as np
import matplotlib as mpl; mpl.use('agg')
import matplotlib.pyplot as plt


def ek_selected(filephmat, filbassel='BASSEL'):
    '''
    Extract energies and k-list that are selected in namd simulation.
    This function will use function read_ephmath5 below.

    Parameters:
    filephmat: string, file name or path of PERTURBO output file.
    filbassel: string, file name or path of BASSEL file of namd simulation.

    Returns: two ndarrays, energies and k-list arrays, in forms of en[nbas]
             and kpts[nbas,3], respectively.
    '''

    bassel  = np.loadtxt(filbassel, dtype=int, skiprows=1) - 1
    en_tot   = read_ephmath5(filephmat, igroup=0, idset=3)
    kpts_tot = read_ephmath5(filephmat, igroup=0, idset=1)

    en = en_tot[bassel[:,0], bassel[:,1]]
    kpts = kpts_tot[bassel[:,0]]

    return en, kpts



def plot_tdshprop(filshps, ptype=1, tdksen=None, figname='tdshp.png'):
    '''
    This function load datas from SHPROP.xxx files,
    and plot average evolution of energy & surface hopping proportion of
    electronic states.

    Parameters:
    filshps: a list of strings, file names if SHPROP.xxx files. such as
             ['SHPROP.1', 'SHPROP.5'].
    lplot  : integer, fig type to plot. 0: do not output figure; 1: plot
             average proportions; 2: plot average energy evolution with
             proportions.
    tdksen : ndarray, time-dependent KS energies in forms of 
             tdksen[ntsteps, nbands]
    figname: string, file name of output figure.

    Returns: ndarray, average data of SHPROP.xxx files.
    '''

    shps = np.array( [ np.loadtxt(filshp) for filshp in filshps ] )
    shp = np.average(shps, axis=0)

    if ptype==0: return shp
        
    figsize_x = 4.8
    figsize_y = 3.2 # in inches
    namdtime = shp[-1,0]

    fig, ax = plt.subplots()
    fig.set_size_inches(figsize_x, figsize_y)
    mpl.rcParams['axes.unicode_minus'] = False

    if ptype==1:
        ylabel = 'SHPROP'
        ax.plot(shp[:,0], shp[:,2:])
        ax.set_ylim(0.0,1.0)
        ax.set_ylabel(ylabel)
    else:
        cmap = 'hot_r'
        dotsize = 20
        ylabel = 'Energy (eV)'
        ntsteps = tdksen.shape[0]
        nbands = tdksen.shape[1]
        vmin = 0; vmax = np.max(shp[:,2:])
        norm = mpl.colors.Normalize(0,vmax)

        T = np.tile(np.arange(ntsteps), nbands).reshape(nbands,ntsteps).T
        ax.scatter(T, tdksen, s=dotsize, c=shp[:,2:], lw=0, 
                   norm=norm, cmap=cmap)
        ax.plot(shp[:,1], 'r', lw=1, label='Average Energy')

        ax.set_ylabel(ylabel)

    ax.set_xlim(0,namdtime)
    ax.set_xlabel('Time (fs)')

    plt.tight_layout()
    plt.savefig(figname, dpi=400)

    return shp


def read_ephmath5(filname, igroup=-1, idset=-1, dset=""):
    '''
    Read informations about e-ph coupling from PERTURBO output h5 file.

    Parameters:
    filname: string, file name or path of PERTURBO output file.
    igroup : integer, which group, start from 0.
    idset  : integer, which data set, start from 0.
    dset   : string, data set name or path, only used when igroup and idset
             are not provided.

    Returns: ndarray, dataset value.
    '''

    if ( igroup > -1  and idset > -1 ):

        group_list = ['el_ph_band_info', 'g_ephmat_total_meV']
        dset_list  = [
            'informations', 'k_list', 'q_list', 'el_band_eV', # 0~3
            'ph_disp_meV', 'phmod_ev_r', 'phmod_ev_i',        # 4~6
            'lattice_vec_angstrom', 'atom_pos', 'mass_a.u.']  # 7~9
        if igroup==0:
            dset_name = dset_list[idset]
        else:
            tag = ['r', 'i']
            dset_name = 'g_ik_' + tag[idset%2] + '_%d'%(idset//2+1)

        f = h5py.File(filname, 'r')
        group = f[group_list[igroup]]
        return group[dset_name].value

    elif (dset!=''):
        path = dset.split('/')[1:]
        f = h5py.File(filname, 'r')
        for item in path:
            dset = f[item]
            f = dset
        return dset.value

    else:
        print("\nNeed input correct args: \'igroup\' & \'idset\' or " \
              "\'dset\'!\n")
        return None


if __name__=='__main__':
    intro = "\nThis is namd postprocessing module, you can add this file to " \
        "your PYTHONPATH,\nand import \"postnamd\" in your python scripts.\n" \
        "\nFor more informations, you can email to zhenfacn@gmail.com.\n"
    print(intro)
