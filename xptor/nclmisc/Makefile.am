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

noinst_LTLIBRARIES = libgacodexptornclmisc.la
libgacodexptornclmisc_la_SOURCES = \
  fluxav_efit.f \
  fluxav_vmec.f \
  geom_circ.f \
  geom_efit.f \
  geom_vmec.f \
  get_pro_xptor.f \
  read_efit_eqdsk.f \
  read_eqdsk.f \
  read_vmec_wout.f \
  store_efit.f \
  track_mod.f

IOINC=@NETCDF_INC@
IOLIBS=@NETCDF_ALIBS@

nclmisc_includes= 
#-I $(top_builddir)/xptor/models \
#  -I $(top_builddir)/tglf/src \
#  -I $(top_builddir)/neo/src 


## To get to work with optimization reducing compiler wrappers
## or simply to get .f to work with new automake
AM_LIBTOOLFLAGS = --tag=FC
AM_FFLAGS=$(FC_FIXED_FLAG) $(FC_DBL_FLAG) $(FC_PIC_FLAG) $(nclmisc_includes) $(IOINC) $(IONLIBS)
AM_FCFLAGS=$(FC_FREE_FLAG) $(FC_DBL_FLAG) $(FC_PIC_FLAG) $(nclmisc_includes) $(IOINC) $(IONLIBS)

CLEANFILES = *.mod *.o *.lo *.lst *.la
