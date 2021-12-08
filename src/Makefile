#-------------------------------------------------------------------------------
# defaults
#-------------------------------------------------------------------------------
FC= ifort -assume byterecl
# FC= gfortran
#FFLAGS= -g -O2
FFLAGS= -g
MAKE = make

#path to HDF5 library
IFLAGS = -I/cluster/apps/hdf5/hdf5-1.12.0/build64/include
LFLAGS = -L/cluster/apps/hdf5/hdf5-1.12.0/build64/lib/ -lhdf5 -lhdf5_fortran
# IFLAGS = -I/public/apps/qe/q-e-qe-6.5_epc/perturbo/h5lib/include
# LFLAGS = -L/public/apps/qe/q-e-qe-6.5_epc/perturbo/h5lib/lib -lhdf5 -lhdf5_fortran
FFLAGS += $(IFLAGS)
# FFLAGS += $(LFLAGS)
SPGLIB += $(LFLAGS)
#-------------------------------------------------------------------------------
# Src
#-------------------------------------------------------------------------------

SRC = prec.f90 lattice.f90 wave.f90 fileio.f90 couplings.f90 epcoup.f90 hamil.f90 \
   	TimeProp.f90 SurfHop.f90  main.f90


OBJ = $(SRC:.f90=.o)
EXE = namd-epc

#-------------------------------------------------------------------------------
# Suffix rules
#-------------------------------------------------------------------------------
.SUFFIXES: .o .f90
.f90.o:
	$(FC) $(FFLAGS) -c $<

#-------------------------------------------------------------------------------
# Targets
#-------------------------------------------------------------------------------
tdm:	$(OBJ)
	$(FC) $(FFLAGS) -o $(EXE) $(OBJ) $(SPGLIB)  

clean:
	rm -f *.mod *.a namd
	rm -f $(OBJ) $(EXE)
tar:
	tar -czvf namd.tgz *.f90 Makefile
tag:
	ctags *.f90
