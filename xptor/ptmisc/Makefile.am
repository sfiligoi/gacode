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

noinst_LTLIBRARIES = libgacodexptorptmisc.la
libgacodexptorptmisc_la_SOURCES = \
  ave_fields.f \
  average7_2d.f \
  calc_stats.f \
  delay.f \
  derivedexpprofiles.f \
  derivedmodlocal.f \
  derivedmodprofiles.f \
  filter_dv.f \
  gridsmooth.f \
  init.f \
  ptorupdate.f \
  setup.f \
  timeprofiles_dv.f \
  tpedscale.f \
  transfer_exp.f \
  update_src.f \
  write_ascii.f \
  write_glf.f \
  write_netcdf.f \
  write_ufile.f

IOINC=@NETCDF_INC@
IOLIBS=@NETCDF_ALIBS@

ptmisc_includes= 
#-I $(top_builddir)/xptor/models \
#  -I $(top_builddir)/tglf/src \
#  -I $(top_builddir)/neo/src 


## To get to work with optimization reducing compiler wrappers
## or simply to get .f to work with new automake
AM_LIBTOOLFLAGS = --tag=FC
AM_FFLAGS=$(FC_FIXED_FLAG) $(FC_DBL_FLAG) $(FC_PIC_FLAG) $(ptlmisc_includes) $(IOINC) $(IONLIBS)
AM_FCFLAGS=$(FC_FREE_FLAG) $(FC_DBL_FLAG) $(FC_PIC_FLAG) $(ptlmisc_includes) $(IOINC) $(IONLIBS)

CLEANFILES = *.mod *.o *.lo *.lst *.la
