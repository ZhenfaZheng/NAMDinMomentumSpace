#!/usr/bin/env python

import math, os, sys
import numpy as np

def main():
    
    nkx = 90
    nky = 90
    nkz = 1

    nks = nkx * nky * nkz

    if (nks <= 1e5):
        nks_p = int( 1e5 / nks) * 100
    elif (nks <= 1e7):
        nks_p = int( 1e7 / nks)
    else:
        print("\nError: Too many k-points (nks > 10^7) !!!\n")
        sys.exit(1)

    nparts = int( math.ceil(float(nks)/nks_p) )
    lenpath = int( np.floor(np.log10(nparts)) ) + 1

    kpoints = kgenerate(nkx, nky, nkz)

    for ip in range(nparts): 

        kst = ip * nks_p
        kend = ( kst + nks_p ) if (ip<nparts-1) else nks
        kpts_p = kpoints[kst:kend,:]
        path = '{:0>{width}d}P'.format(ip+1, width=lenpath)
        if not os.path.exists(path):
            os.makedirs(path)
        filname = path + '/eph.kpt'
        savekpts(filname, kpts_p)
        filname = path + '/eph.qpt'
        savekpts(filname, kpoints)
        # print(path, filname)


def kgenerate(nkx, nky, nkz):

    kstepx = 1.0 / nkx
    kstepy = 1.0 / nky
    kstepz = 1.0 / nkz
    kcenter = np.array([0.0, 0.0, 0.0])
    kpoints = []

    for i in range(nkx):
        offsetx = kstepx * i
        kx = kcenter[0] + offsetx
        for j in range(nky):
            offsety = kstepy * j
            ky = kcenter[1] + offsety
            for k in range(nkz):
                offsetz = kstepz * k
                kz = kcenter[2] + offsetz
                kpoints.append([kx,ky,kz])

    return np.array(kpoints)


def savekpts(filname, kpoints, weight=None):
    
    nks = kpoints.shape[0]

    if weight is None:
        weight = np.ones((nks,1))

    header = '%d'%nks
    data = np.append(kpoints, weight, axis=1)

    np.savetxt(filname, data,
               header ='%d'%nks, comments='',
               fmt='%15.8f%15.8f%15.8f%6d')


if __name__=='__main__':
    main()
