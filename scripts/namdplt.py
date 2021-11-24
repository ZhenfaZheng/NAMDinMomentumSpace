#!/usr/bin/env python

import numpy as np
import postnamd as pn
from glob import glob

def main():

    coup = pn.read_couple(ctype=1)
    coup_av = np.average(np.abs(coup), axis=0)
    pn.plot_couple(coup_av)

    fileptxt = 'EPTXT'
    figname = 'COUPLE_corr.png'
    coup = pn.read_couple(fileptxt, ctype=1)
    coup_av = np.average(np.abs(coup), axis=0)
    pn.plot_couple(coup_av, figname)

    tag = 'SHPROP'
    fileig = 'EIGTXT'
    figname = 'TDEN.png'
    filshps = glob(tag+'.*')
    en = np.loadtxt(fileig)[0,:]
    pn.tdshprop(filshps, lplot=2, ksen=en, figname=figname)


if __name__=='__main__':
    main()
