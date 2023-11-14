#!/usr/bin/env python

import os
import numpy as np
import postnamd as pn


# select initial states with energy
# satisfied | en - inien | < de
inien = 1.0
de = 0.01

Ef = -4.5203
inien += Ef

# select initial states around specified k-point center,
# and satisfies | kpt - kcenter | < dk
kcenter = [0.333, 0.333, 0.0]
dk = 0.3

# select initial time between
# init_time_min ~ init_time_max randomly.
init_time_min = 1
init_time_max = 200

#=====================================================================#

kcenter = np.array(kcenter)
iniemin = inien - de
iniemax = inien + de

# select initial states within energy range.
inp = pn.read_inp('./inp')
en_tot, kpts_tot = pn.read_ektot(inp)
nks = en_tot.shape[0]
nbs = en_tot.shape[1]
index = []
for ik in range(nks):
    kpt = kpts_tot[ik]
    if (np.linalg.norm(kpt-kcenter)<dk):
        for ib in range(nbs):
            en = en_tot[ik, ib]
            if (np.abs(en-inien)<de):
                index.append([ik,ib])
index = np.array(index) + 1

emin = float(inp['EMIN'])
emax = float(inp['EMAX'])
if (iniemin < emin or iniemax > emax):
    print("\nWarning: Initial energy range exceeds EMIN ~ EMAX in inp!\n")

# read number of initial states and determine nsample.
nbas = index.shape[0]
if ('NINIBS' in inp):
    ninibs = int( inp['NINIBS'] )
else:
    ninibs = 1
if (nbas<ninibs):
    print("\nNot enough states between%.2f ~ %.2f eV!\n"%(iniemin, iniemax))
    os._exit(0)
nsample = int( nbas // ninibs )

# generate initial times randomly.
T = np.arange(init_time_min, init_time_max + 1)
np.random.shuffle(T)
niniT = T.shape[0]

if (nsample>niniT):
    nsample = niniT
    print("\nToo many states between %.2f ~ %.2f eV!"%(iniemin, iniemax))
    print("Select first %d states. "%(nsample*ninibs) + \
          "Please set NSAMPLE = %d or less in inp.\n"%niniT)
else:
    print("\nPlease set NSAMPLE = %d or less in inp.\n"%nsample)

# generate data of INICON and save them.
ini = np.empty([nsample,1+ninibs*2])
ini[:,0] = T[:nsample]
ini[:,1:] = index[:nsample*ninibs, :].reshape((nsample,ninibs*2))

np.savetxt('INICON', ini, fmt='%8d')
