#-------------------------------------------------------------------------------
# defaults
#-------------------------------------------------------------------------------
#FC= gfortran
FC= ifort -assume byterecl
FFLAGS= -g -O2
MAKE = make

#path to HDF5 library
IFLAGS = -I/public/apps/qe/hdf5/1.12.1_intel/include
LFLAGS = -L/public/apps/qe/hdf5/1.12.1_intel/lib -lhdf5 -lhdf5_fortran
FFLAGS += $(IFLAGS)
SPGLIB += $(LFLAGS)
#-------------------------------------------------------------------------------
# Src
#-------------------------------------------------------------------------------

SRC = prec.f90 lattice.f90 wave.f90 fileio.f90 couplings.f90 epcoup.f90 hamil.f90 \
   	TimeProp.f90 SurfHop.f90 SurfHopOTF.f90  main.f90


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
	tar -czvf namd-epc.tgz *.f90 Makefile
tag:
	ctags *.f90
