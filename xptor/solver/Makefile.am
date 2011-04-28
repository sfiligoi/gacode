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

noinst_LTLIBRARIES = libgacodexptorsolver.la
libgacodexptorsolver_la_SOURCES = \
  advancefields.f \
  neo_flows.f \
  norm_dv.f \
  perturb_dv.f \
  ptor_dv.f \
  ptorinit_dv.f \
  tglf_startup.f \
  trcoef_dv.f \
  vtrans_dv.f

IOINC=@NETCDF_INC@
IOLIBS=@NETCDF_ALIBS@

solver_includes= \
  -I $(top_builddir)/tglf/src \
  -I $(top_builddir)/neo/src


## To get to work with optimization reducing compiler wrappers
## or simply to get .f to work with new automake
AM_LIBTOOLFLAGS = --tag=FC
AM_FFLAGS=$(FC_FIXED_FLAG) $(FC_DBL_FLAG) $(FC_PIC_FLAG) $(solver_includes) $(IOINC) $(IONLIBS)
AM_FCFLAGS=$(FC_FREE_FLAG) $(FC_DBL_FLAG) $(FC_PIC_FLAG) $(solver_includes) $(IOINC) $(IONLIBS)

CLEANFILES = *.mod *.o *.lo *.lst *.la
