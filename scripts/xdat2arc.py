#!/usr/bin/env python

import ase.io

xdatcar = ase.io.iread(filename='XDATCAR', format='vasp-xdatcar')
ase.io.write('MD.arc',images=list(xdatcar),format='dmol-arc') 
