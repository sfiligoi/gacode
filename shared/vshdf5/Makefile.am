## ##########################################################################
##
## File:		Makefile.am
##
## Author: Mark Fahey
##
## Purpose:  Process this file with automake to produce Makefile
##
## Version:	$Id: Makefile.am 2 2009-09-01 23:34:41Z srinath $
##
## ############################################################################


libdir = $(prefix)/lib

noinst_LTLIBRARIES = libvshdf5.la 

libvshdf5_la_SOURCES = hdf5_fcapi.F90

# To get to work with optimization reducing compiler wrappers
# or simply to get .f to work with new automake
AM_LIBTOOLFLAGS = --tag=FC
AM_FFLAGS=$(FC_FIXED_FLAG) $(FC_DBL_FLAG) $(FC_PIC_FLAG) @HDF5_FC_INC@
AM_FCFLAGS=$(FC_FREE_FLAG) $(FC_DBL_FLAG) $(FC_PIC_FLAG) @HDF5_FC_INC@
CLEANFILES = *.mod *.o *.lo *.lst *.la
