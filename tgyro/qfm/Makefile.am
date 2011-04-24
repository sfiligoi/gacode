## ##########################################################################
##
## File:		Makefile.am
##
## Author: Mark Fahey
##
## Purpose:  Process this file with automake to produce Makefile
##
## Version:	$Id: Makefile.am 4 2009-09-02 00:04:18Z srinath $
##
## ############################################################################


libdir = $(prefix)/lib
noinst_LTLIBRARIES = libgacodetgyroqfm.la
libgacodetgyroqfm_la_SOURCES = \
  difftd.f \
  disp9t.f \
  qfm_interface.f90 \
  qfm_run.f90 \
  qfm_sub.f \
  r8tomsqz.f

# To get to work with optimization reducing compiler wrappers
# or simply for directories with mixed F77 and FC
AM_LIBTOOLFLAGS = --tag=FC
AM_FFLAGS=$(FC_FIXED_FLAG) $(FC_DBL_FLAG)
AM_FCFLAGS=$(FC_FREE_FLAG) $(FC_DBL_FLAG) $(FC_PIC_FLAG)
CLEANFILES = *.mod *.o *.lo *.lst *.la
