#!/usr/bin/env python

import numpy as np
import postnamd as pn
from glob import glob

def main():

    Eref = -4.52862

    coup = pn.read_couple(ctype=1)
    coup_av = np.average(np.abs(coup), axis=0)
    pn.plot_couple(coup_av)

    fileptxt = 'EPTXT'
    figname = 'COUPLE_sh.png'
    coup = pn.read_couple(fileptxt, ctype=2)
    coup_av = np.average(np.abs(np.sum(coup, axis=0)), axis=0)
    pn.plot_couple(coup_av, figname)

    tag = 'SHPROP'
    fileig = 'EIGTXT'
    figname = 'TDEN.png'
    filshps = glob(tag+'.*')
    filephmat = '../graphene_ephmat_p1.h5'
    shp = pn.readshp(filshps)
    en, kpts = pn.ek_selected(filephmat)
    # Make sure shp & ksen have same Eref!!!
    # pn.plot_namd_3D(shp, kpts[:,:2], en, Eref=Eref)
    pn.plot_tdprop(shp, Eref, lplot=2, ksen=en, figname=figname)
    en_tot = pn.read_ephmath5(filephmat, igroup=0, idset=3)
    kpts_tot = pn.read_ephmath5(filephmat, igroup=0, idset=1)
    pn.plot_namd_3D(shp, kpts[:,:2], en, kpts_tot[:,:2], en_tot[:,1], Eref)

    tag = 'PSICT'
    fileig = 'EIGTXT'
    figname = 'TDPSI.png'
    filshps = glob(tag+'.*')
    shp = pn.readshp(filshps)
    pn.plot_tdprop(shp, Eref, lplot=2, ksen=en, figname=figname)


if __name__=='__main__':
    main()
