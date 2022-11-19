#!/usr/bin/env python

import os
import numpy as np
import postnamd as pn


# select initial states
# between iniemin ~ iniemax
iniemin = 0.59
iniemax = 0.61
Ef = -4.52862
iniemin += Ef
iniemax += Ef

# select initial time between
# init_time_min ~ init_time_max randomly.
init_time_min = 1
init_time_max = 200

#=====================================================================#

inp = pn.read_inp('./inp')
nparts = int(inp['NPARTS'])
prefix = inp['EPMPREF']
epmdir = inp['EPMDIR']

# read band energies from h5 files.
for ip in range(nparts):
    filepm = os.path.join(epmdir, prefix + '_ephmat_p%d.h5'%(ip+1))
    en_p = pn.read_ephmath5(filepm, igroup=0, idset=3)
    if ip==0:
        en_tot = en_p
    else:
        en_tot = np.vstack((en_tot, en_p))

# select initial states within energy range.
index = np.argwhere((en_tot>iniemin) & (en_tot<iniemax)) + 1

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
