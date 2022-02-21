#!/usr/bin/env python

import numpy as np
import postnamd as pn
from glob import glob


def main():

    Eref = -4.52862
    filephmat = './graphene_ephmat_p1.h5'
    en_tot = pn.read_ephmath5(filephmat, igroup=0, idset=3)
    kpts_tot = pn.read_ephmath5(filephmat, igroup=0, idset=1)
    kpath = np.array([
        [0.33333, 0.33333, 0.00000],
        [0.00000, 0.00000, 0.00000],
        [0.50000, 0.00000, 0.00000],
        [0.33333, 0.33333, 0.00000]
        ])
    bvec = np.array([
        [1.00000, 0.00000, 0.00000],
        [0.00000, 1.00000, 0.00000],
        [0.00000, 0.00000, 1.00000]
        ])

    k_tag = select_kpts_on_path(kpts_tot, kpath)
    # np.savetxt('kpts.dat', kpts_tot, fmt='%8.4f')
    # np.savetxt('ktag.dat', k_tag, fmt='%5d')
    loc_on_kpath(kpts_tot, k_tag, kpath, bvec)


def loc_on_kpath(kpts, k_tag, kpath, bvec):
    
    npath = kpath.shape[0] - 1
    segment_lenth = np.zeros(npath)
    for ipath in range(npath):
        dk = kpath[ipath+1] - kpath[ipath]
        dk_cart = np.zeros(3)
        for iax in range(3):
            dk_cart += dk[iax]*bvec[iax,:]
        segment_lenth[ipath] = np.linalg.norm(dk_cart)
    print(segment_lenth)

def select_kpts_on_path(kpts, kpath):

    nk = kpts.shape[0]
    npath = kpath.shape[0] - 1

    if npath==0:
        print("\nkpath must contain 2 k-points at least!\n")
        return

    pathes = np.zeros( (npath, 2, 3) )
    pathes[:,0,:] = kpath[:npath, :]
    pathes[:,1,:] = kpath[1:, :]
    dks = pathes[:,1,:] - pathes[:,0,:]

    #k_tag = np.zeros((nk, npath), dtype='int32')
    k_tag = []
    for ik in range(nk):
        kpt = kpts[ik,:]
        for ipath in range(npath):
            zz = k_on_path(kpt, pathes[ipath], dks[ipath])
            if (zz):
                #k_tag[ik, ipath] = 1
                k_tag.append([ik+1,ipath+1])

    k_tag = np.array(k_tag, dtype='int32')
    return k_tag


def k_on_path(kpt, path, dk, norm=1.0e-4):

    # whether kpt in the range of path
    for iax in range(3):
        if ( dk[iax] > 0 ):
            if ( kpt[iax]<(path[0,iax]-norm) or \
                 kpt[iax]>(path[1,iax]+norm) ):
                return False
        else:
            if ( kpt[iax]>(path[0,iax]+norm) or \
                 kpt[iax]<(path[1,iax]-norm) ):
                return False

    # (x-x1) * dy .vs. (y-y1) * dx
    a = (kpt[0]-path[0,0]) * dk[1]
    b = (kpt[1]-path[0,1]) * dk[0]
    if ( abs(b-a)>norm ): return False

    # (y-y1) * dz .vs. (z-z1) * dy
    a = (kpt[1]-path[0,1]) * dk[2]
    b = (kpt[2]-path[0,2]) * dk[1]
    if ( abs(b-a)>norm ): return False

    # (x-x1) * dz .vs. (z-z1) * dx
    a = (kpt[0]-path[0,0]) * dk[2]
    b = (kpt[2]-path[0,2]) * dk[0]
    if ( abs(b-a)>norm ): return False

    # Not including end point of path
    if (np.sum(np.abs(kpt-path[1,:]))<norm): return False

    return True


if __name__=='__main__':
    main()
