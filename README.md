
# NAMD in Momentum Space (Hefei-NAMD-EPC)

[![arXiv shield](https://img.shields.io/badge/arXiv-2210.00529-red.svg?style=flat)](https://doi.org/10.38550/arXiv.2210.00529)

## Contents
- [Overview](#overview)
- [Installation Guide](#installation-guide)
- [Tutorials of Runing Hefei-NAMD-EPC](#tutorials-of-runing-hefei-namd-epc)
- [Citation](#citation)


# Overview

This is a program to perform non-adiabatic molecular dynamics (NAMD) simultion
in momentum space, where the non-adiabatic coupling (NAC) are replaced by
electron-phonon (*e-ph*) coupling (EPC).

The codes are based on [Hefei-NAMD](https://github.com/QijingZheng/Hefei-NAMD)
developed by [Qijing Zheng](http://staff.ustc.edu.cn/~zqj), et al.


# Installation Guide

The user can get the source files of this program from Github
```
git clone https://github.com/ZhenfaZheng/NAMDinMomentumSpace.git
```
Then, go into folder of source files
```
cd src
```
To compile this program, you should specify the path to your HDF5 library.
Open the "Makefile"
```
vim Makefile
```
Modify the two flags "IFLAGS" and "LFLAGS"
```
IFLAGS = -I/\${path-to-your-hdf5-dir}/include
LFLAGS = -I/\${path-to-your-hdf5-dir}/include -lhdf5 -lhdf5_fortran
```
After modifying the "Makefile", compile the Hefei-NAMD-EPC
```
make clean
make
```
You will get an executable file "namd-epc" for simulations of NAMD in momentum
space.


# Tutorials of Runing Hefei-NAMD-EPC

1. Calculate *e-ph* matrix elements using PERTURBO

Befor carrying out NAMD simulations in momentum space, the user needs to
calculate *e-ph* matrix elements $g\_{mn\nu}(\mathbf{k}, \mathbf{q})$ using a
modified version PERTURBO. Therefore, first of all, the user needs to install
and learn using PERTURBO. We recommend users to visit the home page of PERTURBO
([https://perturbo-code.github.io](https://perturbo-code.github.io)) to get the
source code and tutorials of PERTURBO 1.0 which compatible with the 6.5 and 6.3
version Quantum Espresso (QE).

Then, replace the source file "calc\_ephmat.f90" in path of "pert-src/" with
our modified file in the path of "patch" before compiling PERTURBO. It should
be noted that this patch is modified by Prof. Jin-Jian Zhou and us together,
and it will output files containing *e-ph* information (named
"\${prefix}\_ephmat\_p1.h5", "\${prefix}\_ephmat\_p2.h5",
"\${prefix}\_ephmat\_p2.h5", ...) additionally when run "perturbo.x" with
"calc\_mode = 'ephmat'".

2. Prepare input files

To run Hefei-NAMD-EPC, the user should prepare 3 kinds of files:
"\${prefix}\_ephmat\_p\${ipart}.h5" files containing *e-ph* information from
PERTURBO calculations, "inp" file for input parameters setting, and "INICON"
file for initial conditions setting.

First, copy or link "\${prefix}\_ephmat\_p\${ipart}.h5" files from PERTURBO
calculations
```
mkdir epfiles
cd epfiles
for ipart in {1..8}  # Here, 8 is the No. of parts as an example.
do
ln -sf ${path-to-your-perturbo-work-dir}/${prefix}_ephmat_p${ipart}.h5 .
done
cd ..
```

Then, generate "inp" file, below is an example
```
&NAMDPARA
  BMIN       = 1
  BMAX       = 2
  KMIN       = 1
  KMAX       = 40
  EMIN       = -4.6
  EMAX       = -1.5
  NBANDS     = 2
  NKPOINTS   = 81
  NINIBS     = 1
  Np         = 81

  NSW        = 1200
  POTIM      = 1.0
  TEMP       = 300

  NSAMPLE    = 6
  NAMDTIME   = 1000
  NELM       = 1000
  NTRAJ      = 10000
  LHOLE      = .F.
  LCPEXT     = .F.

  LEPC       = .T.
  LBASSEL    = .F.
  LARGEBS    = .T.
  LSORT      = .T.
  EPCTYPE    = 1
  NPARTS     = 1
  EPMDIR     = './'
  EPMPREF    = 'graphene'
/
```

The last file we need is "INICON" file, for example
```
   11     4     2
   23    10     1
   47    17     2
   51    21     2
   85    34     2
   93    38     2
```
Each line of "INICON" represents initial condition of each sample, the first
column indicates initial time, the second and third column indicate
$\mathbf{k}$ and band indices of initial electronic state.


3. Run Hefei-NAMD-EPC

Copy executable file "namd-epc" to your work diretory, then execute
```
./namd-epc
```

4. Output files

After Hefei-NAMD-EPC calculation, we obtain files as bellow:

BASSEL (Belected basises in NAMD simulation)

EIGTXT (Eigen energy of each basis)

EPTXT, EPECTXT, EPPHTXT (*e-ph* couplings)

PHPROP.\${imode} (Evolution of phonon number for each mode)

PSICT.\${initime} (Evolution of expanding coefficients for each sample)

SHPROP.\${initime} (Evolution of population for each sample)

5. Data processing

We provide python scripts to do the data processing, just copy the
"postnamd.py" and "namdplt.py" from the folder "scripts" to your work diretory,
and execute
```
python namdplt.py
```
You will obtain figures of results, including: COUPLE.png, COUPLE\_EL.png,
COUPLE\_PH.png, TDEN.png, TDBAND.png, TDPH.png, TDPHEN.png.


# Citation

For theoretical details on the code, we refer the users to the manuscript
accompying the source code:

> Z. Zheng, Y. Shi, J.J. Zhou, O. V. Prezhdo, Q. Zheng, J. Zhao. *"Ab initio
real-time quantum dynamics of charge carriers in momentum space"*.
[arXiv:2210.00529 (2022)](https://doi.org/10.38550/arXiv.2210.00529)

When using results from this program in your publications, please cite the
paper given above and acknowledge the use of the code.

It would also be appropriate to cite the original articles:

> Q. Zheng, W. A. Saidi, Y. Xie, Z. Lan, O. V. Prezhdo, H. Petek, J. Zhao.
"*Phonon-Assisted Ultrafast Charge Transfer at van der Waals Heterostructure
Interface*".
[Nano Lett 17, 6435-6442, (2017)](https://doi.org/10.1021/acs.nanolett.7b03429)

> Q. Zheng, W. Chu, C. Zhao, L. Zhang, H. Guo, Y. Wang, X. Jiang, J. Zhao.
"*Ab initio nonadiabatic molecular dynamics investigations on the excited
carriers in condensed matter systems*". 
[WIREs Comput. Mol. Sci. 9:e1411, (2019)](https://doi.org/10.1002/wcms.1411)


