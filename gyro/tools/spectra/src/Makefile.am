## ############################################################
##
## File:  Makefile.am
##
## Purpose:  
##
## $Id: Makefile.am 38 2010-01-11 21:51:29Z srinath $
## ############################################################


# To get to work with optimization reducing compiler wrappers
#IOINC=@HDF5_INC@ @FCIOWRAPPERS_INC@ @NETCDF_INC@
#IOINC=@NETCDF_INC@
IOINC=
#IOLIBS=@FCIOWRAPPERS_LIBS@ @HDF5_LIBS@ @NETCDF_LIBS@
IOLIBS=@NETCDF_LIBS@
AM_LIBTOOLFLAGS = --tag=FC
AM_FCFLAGS=$(FC_FREE_FLAG) $(FC_DBL_FLAG) $(FC_PIC_FLAG) $(IOINC) -I.
AM_FFLAGS=$(FC_FIXED_FLAG) $(FC_DBL_FLAG) $(FC_PIC_FLAG) $(IOINC) -I.


bin_PROGRAMS = spec_xy spec_n spec_i
spec_xy_SOURCES = \
  input_data.f90 \
  read_input.f90 \
  per_deriv.f90 \
  spec_xy.f90

spec_n_SOURCES = \
  input_data.f90 \
  read_input.f90 \
  per_deriv.f90 \
  spec_n.f90

spec_i_SOURCES = \
  input_data.f90 \
  read_input.f90 \
  per_deriv.f90 \
  spec_i.f90
   

MAINTAINERCLEANFILES    = Makefile.in  *.mod *.o

   
