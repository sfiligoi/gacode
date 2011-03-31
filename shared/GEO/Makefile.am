## ##########################################################################
##
## File:		Makefile.am
##
## Author: 
##
## Purpose:  Process this file with automake to produce Makefile
##
## Version:	$Id: Makefile.am 2 2009-09-01 23:34:41Z srinath $
##
## ############################################################################


libdir = $(prefix)/lib
noinst_LTLIBRARIES = libgacodegeo.la
libgacodegeo_la_SOURCES = \
  GEO_interface.f90 \
  GEO_alloc.f90 \
  GEO_do.f90 \
  GEO_interp.f90 \
  GEO_write.f90

# To get to work with optimization reducing compiler wrappers
AM_LIBTOOLFLAGS = --tag=FC
AM_FFLAGS=$(FC_FIXED_FLAG) $(FC_DBL_FLAG) $(FC_PIC_FLAG)
AM_FCFLAGS=$(FC_FREE_FLAG) $(FC_DBL_FLAG) $(FC_PIC_FLAG)
CLEANFILES = *.mod *.o *.lo *.lst *.la
