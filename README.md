
# NAMD in Momentum Space (Hefei-NAMD-EPC)

[![arXiv shield](https://img.shields.io/badge/arXiv-2210.00529-red.svg?style=flat)]
(https://doi.org/10.48550/arXiv.2210.00529)

# Overview

This is a program to perform non-adiabatic molecular dynamics (NAMD) simultion
in momentum space, where the non-adiabatic coupling (NAC) are replaced by
electron-phonon (*e-ph*) coupling (EPC).

The codes are based on [Hefei-NAMD](https://github.com/QijingZheng/Hefei-NAMD)
developed by [Qijing Zheng](http://staff.ustc.edu.cn/~zqj), et al.

# Citation

For more details on the code, we refer the users to the manuscript accompying
the source code:

> Z. Zheng, Y. Shi, J.J. Zhou, O. V. Prezhdo, Q. Zheng, J. Zhao, *"Ab initio
real-time quantum dynamics of charge carriers in momentum space"*,
[arXiv:2210.00529 (2022)](https://doi.org/10.48550/arXiv.2210.00529)

When using results from this program in your publications, please cite the
paper given above and acknowledge the use of the code.

It would also be appropriate to cite the original articles:


# Tutorials of runing Hefei-NAMD-EPC

## Calculate *e-ph* matrix elements using PERTURBO

Befor carrying out NAMD simulations in momentum space, the user needs to
calculate *e-ph* matrix elements $g\_{mn\nu}(\mathbf{k}, \mathbf{q})$ using a
modified version PERTURBO. Therefore, first of all, the user needs to install
and learn using PERTURBO. We recommend users to visit the home page of PERTURBO
[https://perturbo-code.github.io](https://perturbo-code.github.io) to get the
source code and tutorials of PERTURBO 1.0 which compatible with the 6.5 and 6.4
version Quantum Espresso (QE).

Then, replace the source file "calc\_ephmat.f90" in path of "pert-src/" with
our modified file in the path of "patch" before compiling PERTURBO. It should
be noted that this patch is modified by Prof. Jin-Jian Zhou and us together,
and it will output files containing *e-ph* information (named
"\\${prefix}\_ephmat\_p1.h5", "\${prefix}\_ephmat\_p2.h5",
"\${prefix}\_ephmat\_p2.h5", ...) additionally when run "perturbo.x" with
"calc\_mode = 'ephmat'".
