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

noinst_LTLIBRARIES = libgacodexptor_SRC.la
libgacodexptor_SRC_la_SOURCES = \
  neucg.f \
  neuden.f \
  neut.f \
  palpha_dv.f \
  pohmic_dv.f \
  prad_dv.f \
  sigvth.f

IOINC=@NETCDF_INC@
IOLIBS=@NETCDF_ALIBS@

source_includes= \
  -I $(top_builddir)/xptor/src \
  -I $(top_builddir)/tglf/src \
  -I $(top_builddir)/neo/src


## To get to work with optimization reducing compiler wrappers
## or simply to get .f to work with new automake
AM_LIBTOOLFLAGS = --tag=FC
AM_FFLAGS=$(FC_FIXED_FLAG) $(FC_DBL_FLAG) $(FC_PIC_FLAG) $(source_includes) $(IOINC) $(IONLIBS)
AM_FCFLAGS=$(FC_FREE_FLAG) $(FC_DBL_FLAG) $(FC_PIC_FLAG) $(source_includes) $(IOINC) $(IONLIBS)

CLEANFILES = *.mod *.o *.lo *.lst *.la
