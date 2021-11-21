#!/usr/bin/env python3

import numpy as np
import matplotlib as mpl; mpl.use('agg')
import matplotlib.pyplot as plt

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


if __name__=='__main__':
    intro = "\nThis is namd postprocessing module, you can add this file to " \
        "your PYTHONPATH,\nand import \"postnamd\" in your python scripts.\n" \
        "\nFor more informations, you can email to zhenfacn@gmail.com.\n"
    print(intro)
