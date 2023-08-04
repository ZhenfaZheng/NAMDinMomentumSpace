#!/usr/bin/env python

import numpy as np
import os, math, h5py
from scipy.optimize import curve_fit


def read_inp(infile='inp'):
    '''
    Read input parameters.

    Parameters:
    infile: string, coupling file.

    Returns: dictionary, the keys and values are all strings.
    '''

    text = [line for line in open(infile) if line.strip()]

    inp = {}
    for line in text:
        if (line[0]=='&' or line[0]=='/' or line[0]=='!'):
            continue
        temp = line.split('!')[0].split('=')
        key = temp[0].strip()
        value = temp[1].strip().strip('\'').strip('\"')
        inp[key] = value

    return inp


def read_couple(filcoup='NATXT', inp=None, ctype=0):
    '''
    This function loads data from NATXT file.

    Parameters:
    filcoup: string, coupling file.
    inp    : dictionary, input parameters.
    ctype  : integer, different forms of NATXT files.
             0: origin type, NAC data are real number;
             1: NAC are complex, and restore in two real numbers;
             2: EPC divide into NM (this number stored in first
             line) parts, each part has same form with type 1;

    Returns: ndarray, coupling data in forms of coup[nsw-1, nb, nb]
             or coup[NM, nsw-1, nb, nb] (type 2)
    '''

    if (inp is not None):
        largeBS = inp['LARGEBS'].strip('.')[0]
        if (largeBS=='T' or largeBS=='t'):
            ctype = 0
        else:
            ctype = 1

    if ctype==0:
        coup = np.loadtxt(filcoup)
        try:
            nt = int(coup.shape[0])
            nb = int( np.sqrt(coup.shape[1]) )
        except IndexError:
            nt = 1
            nb = int( np.sqrt(coup.shape[0]) )

        coup.resize(nt,nb,nb)

    elif ctype==1:
        data = np.loadtxt(filcoup)
        try:
            nt = int(data.shape[0])
            nb = int( np.sqrt(data.shape[1]/2) )
            coup = data[:,0::2] + data[:,1::2]*(1.0j)
        except IndexError:
            nt = 1
            nb = int( np.sqrt(data.shape[0]/2) )
            coup = data[0::2] + data[1::2]*(1.0j)

        coup.resize(nt,nb,nb)

    elif ctype==2:
        data = np.loadtxt(filcoup, skiprows=1)
        with open(filcoup) as f:
            first_line = f.readline()
        nm = int( first_line.split()[0] )
        nb = int( np.sqrt(data.shape[1]/2) )
        nt = int(data.shape[0] / nm)

        coup = data[:,0::2] + data[:,1::2]*(1.0j)
        coup.resize(nm,nt,nb,nb)

    return coup

def read_ektot(inp):
    '''
    Extract total energies and k-list.

    Parameters:
    inp      : dictionary, input parameters.

    Returns: two ndarrays, energies and k-list arrays, in forms of
             en_tot[nks, nb] and kpts_tot[nks,3], respectively.
    '''

    nparts = int(inp['NPARTS'])
    prefix = inp['EPMPREF']
    epmdir = inp['EPMDIR']

    for ip in range(nparts):
        filepm = prefix + '_ephmat_p%d.h5'%(ip+1)
        path = os.path.join(epmdir, filepm)
        en_p   = read_ephmath5(path, igroup=0, idset=3)
        kpts_p = read_ephmath5(path, igroup=0, idset=1)
        if ip==0:
            en_tot = en_p
            kpts_tot = kpts_p
        else:
            en_tot = np.vstack((en_tot, en_p))
            kpts_tot = np.vstack((kpts_tot, kpts_p))

    return en_tot, kpts_tot


def ek_selected(filephmat='', filbassel='BASSEL', inp=None):
    '''
    Extract energies and k-list that are selected in namd simulation.
    This function will use function read_ephmath5 below.

    Parameters:
    filephmat: string, file name or path of PERTURBO output file.
    filbassel: string, file name or path of BASSEL file of namd simulation.
    inp      : dictionary, input parameters.

    Returns: two ndarrays, energies and k-list arrays, in forms of en[nbas]
             and kpts[nbas,3], respectively.
    '''

    if inp is None:
        en_tot   = read_ephmath5(filephmat, igroup=0, idset=3)
        kpts_tot = read_ephmath5(filephmat, igroup=0, idset=1)
    else:
        en_tot, kpts_tot = read_ektot(inp)

    bassel  = np.loadtxt(filbassel, skiprows=1)
    kb_index = np.array(bassel[:, :2], dtype=int) - 1
    en = en_tot[kb_index[:,0], kb_index[:,1]]
    kpts = kpts_tot[kb_index[:,0]]

    return en, kpts

def read_inicon(filinicon='INICON'):
    '''
    This function loads data from INICON file.
    '''
    inicon = np.loadtxt(filinicon)
    return inicon


def readshp(filshps):
    '''
    This function loads data from SHPROP.xxx files.

    Parameters:
    filshps: a list of strings, file names if SHPROP.xxx files. such as
             ['SHPROP.1', 'SHPROP.5']

    Returns: ndarray, average data of SHPROP.xxx files, in forms of
             shp[ntsteps, nb+2].
    '''
    shps = np.array( [ np.loadtxt(filshp) for filshp in filshps ] )
    shp = np.average(shps, axis=0)
    return shp


def readphp(filphps):
    '''
    This function loads data from SHPROP.xxx files.

    Parameters:
    filphps: a list of strings, file names if PHPROP.xxx files. such as
             ['PHPROP.1', 'PHPROP.2']

    Returns: ndarray, total data of PHPROP.xxx files, in forms of
             php[nmodes, ntsteps, nb+2].
    '''
    php = np.array( [ np.loadtxt(filphp) for filphp in filphps ] )
    return php


def read_ephmath5(filname, igroup=-1, idset=-1, dset=""):
    '''
    Read informations about e-ph coupling from PERTURBO output h5 file.

    Parameters:
    filname: string, file name or path of PERTURBO output file.
    igroup : integer, which group, start from 0.
    idset  : integer, which data set, start from 0.
    dset   : string, data set name or path, only used when igroup and idset
             are not provided. e.g., "/el_ph_band_info/k_list".

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
        # return group[dset_name].value
        return group[dset_name][:]

    elif (dset!=''):
        path = dset.split('/')[1:]
        f = h5py.File(filname, 'r')
        for item in path:
            dset = f[item]
            f = dset
        # return dset.value
        return dset[:]

    else:
        print("\nNeed input correct args: \'igroup\' & \'idset\' or " \
              "\'dset\'!\n")
        return None


def calc_rec_vec(a1, a2, a3):
    '''
    Calculate reciprocal lattice vectors.
    Can also calculate basic lattice vectors from reciprocal vectors.

    a1, a2, a3: ndarray, with shape of (3,), basic lattice vectors.

    Returns:
    b1, b2, b3: ndarray, with shape of (3,), reciprocal lattice vectors.
    '''
    b1 = 2 * np.pi * np.cross(a2, a3) / np.dot(a1, np.cross(a2, a3))
    b2 = 2 * np.pi * np.cross(a3, a1) / np.dot(a1, np.cross(a2, a3))
    b3 = 2 * np.pi * np.cross(a1, a2) / np.dot(a1, np.cross(a2, a3))
    return b1, b2, b3


def frac2cart(X, a1, a2, a3):
    '''
    Convert from fractional coordinates to Cartisian coordinates.

    X: ndarray, with shape of (n,3) or (3,), fractional coordinates.
    a1, a2, a3: ndarray, with shape of (3,), basic lattice vectors.

    Returns:
    Y: ndarray, with same shape of X, Cartisian coordinates.
    '''

    A = np.array([a1, a2, a3])
    Y = np.matmul(X, A)
    return Y


def cart2frac(X, a1, a2, a3):
    '''
    Convert from Cartisian coordinates to fractional coordinates.

    X: ndarray, with shape of (n,3) or (3,), Cartisian coordinates.
    a1, a2, a3: ndarray, with shape of (3,), basic lattice vectors.

    Returns:
    Y: ndarray, with same shape of X, fractional coordinates.
    '''

    b1, b2, b3 = calc_rec_vec(a1, a2, a3)
    B = np.array([b1, b2, b3]).T / (2*np.pi)
    Y = np.matmul(X, B)
    return Y


def loc_on_kpath(kpts, k_index, path_index, kpath):
    '''
    Calculate locations of k points on path.

    Parameters:
    kpts      : array, with shape of (nks,3) or (3), k point coordinates.
    k_index   : array, with shape of (nks_onpath), indexes of k points on path.
    path_index: array, with shape of (nks_onpath), indexes of k path segments.
    kpath     : array, with shape of (npath+1, 3), coordinates of begin and
                end points of k path.

    Returns:
    loc: array, with shape of (nks_onpath), locations of k points on path.
    '''
    
    if (len(kpts.shape)==1):
        kpts = np.array([kpts])

    npath = kpath.shape[0] - 1
    dks = np.sum( np.abs(np.diff(kpath, axis=0)), axis=1)
    nks_onpath = k_index.shape[0]

    kpath_loc = np.zeros(npath+1)
    segment_length = np.zeros(npath)
    for ipath in range(npath):
        dk = kpath[ipath+1] - kpath[ipath]
        segment_length[ipath] = np.linalg.norm(dk)
        kpath_loc[ipath+1] = kpath_loc[ipath] + segment_length[ipath]

    loc = np.zeros(nks_onpath)
    for ii in range(nks_onpath):
        ik = k_index[ii]
        ipath = path_index[ii]
        scale = np.sum( np.abs(kpts[ik] - kpath[ipath]) ) / dks[ipath]
        loc[ii] = kpath_loc[ipath] + segment_length[ipath] * scale

    return loc, kpath_loc


def select_kpts_on_path(kpts, kpath, norm=1e-4):
    '''
    Select k points on k path from kpts.

    Parameters:
    kpts : array, with shape of (nks,3) or (3), k point coordinates.
    kpath: array, with shape of (npath+1, 3), coordinates of begin and end
           points of k path.
    norm : float, stard value to determine whether k point is on path.

    Returns:
    k_index   : array, with shape of (nks_onpath), indexes of k points on path.
    path_index: array, with shape of (nks_onpath), indexes of k path segments.
    '''

    if (len(kpts.shape)==1):
        kpts = np.array([kpts])
        nks = 1
    else:
        nks = kpts.shape[0]

    npath = kpath.shape[0] - 1

    if npath==0:
        print("\nkpath must contain 2 k-points at least!\n")
        return

    pathes = np.zeros( (npath, 2, 3) )
    pathes[:,0,:] = kpath[:npath, :]
    pathes[:,1,:] = kpath[1:, :]
    dks = pathes[:,1,:] - pathes[:,0,:]

    #k_tag = np.zeros((nks, npath), dtype='int32')
    k_index = []; path_index = []
    for ik in range(nks):
        kpt = kpts[ik,:]
        for ipath in range(npath):
            lendkpt = False if (ipath<npath-1) else True
            zz = k_on_path(kpt, pathes[ipath], dks[ipath], norm, lendkpt)
            if (zz):
                #k_tag[ik, ipath] = 1
                k_index.append(ik); path_index.append(ipath)

    k_index = np.array(k_index, dtype='int32')
    path_index = np.array(path_index, dtype='int32')
    return k_index, path_index


def k_on_path(kpt, path, dk, norm, lendkpt=False):
    '''
    Determine whether k point is on path.

    Parameters:
    kpt : array, with shape of (3), k point coordinates.
    path: array, with shape of (2, 3), coordinates of begin and end points
          of path.
    dk  : array, with shape of (3), difference between begin and end points
          of path.
    norm: float, stard value to determine whether k point is on path.
    lendkpt: bool, whether to include end point of path

    Returns: True or False.
    '''

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

    if not lendkpt:
        # do not include end point of path
        if (np.sum(np.abs(kpt-path[1,:]))<norm): return False

    return True


def get_Enk(kpath, B, inp):
    '''
    Extract eigen values on kpath.
    '''
    en, kpts= read_ektot(inp)

    b1, b2, b3 = (B[0], B[1], B[2])
    kpts_cart = frac2cart(kpts, b1, b2, b3)
    kpath_cart = frac2cart(kpath, b1, b2, b3)

    k_index, kp_index = select_kpts_on_path(kpts, kpath)
    k_loc, kp_loc = loc_on_kpath(kpts_cart, k_index, kp_index, kpath_cart)

    Enk = np.hstack(( np.array([k_loc]).T, en[k_index,:] ))
    sort = np.argsort(Enk[:,0])

    return Enk[sort,:]


def fit_decaytime(time, ft, ftype=1):
    '''
    Fit decay time.

    Parameters:
    time : ndarray, evolution time, in forms of time[nt].
    ft   : ndarray, time evolution of a specified quantity, in forms of
           ft[nt].
    ftype: integer, fitting type. 1: exponential; 2: gaussian.

    Returns:
    ffit: ndarray, fitting function, in forms of ffit[nt].
    dct : float, decay time.
    '''

    if ftype==1:
        popt, pcov = curve_fit(func_exp, time, ft)
        ffit = func_exp(time, *popt)
        dct = popt[1]
    elif ftype==2:
        popt, pcov = curve_fit(func_gauss, time, ft)
        ffit = func_gauss(time, *popt)
        dct = popt[2]

    return ffit, dct


def func_exp(x, a, b, c):
    '''
    Exponential function.

    Parameters:
    x : ndarray, function variable.
    a, b, c : float, function patameters.
    '''
    return a * np.exp(-x/b) + c


def func_gauss(x, a, b, sigma, c):
    '''
    Gaussian function.

    Parameters:
    x : ndarray, function variable.
    a, b, sigma, c : float, function patameters.
    '''
    return a * np.exp( -0.5 * ((x-b)/sigma)**2 ) + c

def func_fd(en, T):
    '''
    Fermi-dirac distribution

    en: ndarray, function variable, should be subtracted by fermi energy.
    T : float, temperature in unit of K.
    '''
    kbT = 8.61733 * 1.0e-5 * T
    f = 1.0 / ( np.exp(en/kbT) + 1.0 )

    return f


if __name__=='__main__':
    intro = "\nThis is a namd postprocessing module, you can add this file to" \
        " your PYTHONPATH,\nand import \"postnamd\" in your python scripts.\n" \
        "\nFor more informations, you can email to zhenfacn@gmail.com.\n"
    print(intro)
