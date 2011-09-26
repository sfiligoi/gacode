## ############################################################
##
## File:  Makefile.am
##
## Purpose:  
##
## $Id: Makefile.am 38 2010-01-11 21:51:29Z srinath $
## ############################################################


# To get to work with optimization reducing compiler wrappers
IOINC=@NETCDF_INC@ @HDF5_INC@ 
IOLIBS=@NETCDF_ALIBS@ @HDF5_ALIBS@ @Z_LIBS@
AM_LIBTOOLFLAGS = --tag=FC
AM_FCFLAGS=$(FC_FREE_FLAG) $(FC_DBL_FLAG) $(FC_PIC_FLAG) $(IOINC) -I.
AM_FFLAGS=$(FC_FIXED_FLAG) $(FC_DBL_FLAG) $(FC_PIC_FLAG) $(IOINC) -I.


bin_PROGRAMS = gyro2gkv_exec
gyro2gkv_exec_SOURCES = \
  gyro2gkv_globals.f90 \
  read_all.f90 \
  gyro2gkv_exec.f90    

gyro2gkv_exec_LDADD = $(IOLIBS)


EXTRA_DIST = \
  README \
  gread.pro

MAINTAINERCLEANFILES    = Makefile.in  *.mod *.o

   
