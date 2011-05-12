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

noinst_LTLIBRARIES = libgacodexptormodels.la
libgacodexptormodels_la_SOURCES = \
  cgg_solver.f \
  etaw17a.f \
  f77_dcomplx.h \
  forcebal.f \
  glf2d.F \
  glf2d_dv.f \
  glf2d_max.f \
  ifsppl_dv.f \
  ip_chi1.f \
  ip_chi2.f \
  mmm95.f \
  model.f \
  nclass_mod.f \
  neoclass.f \
  neoclass_one.f \
  weiland14.f \
  xptor_neo_map.f 

IOINC=@NETCDF_INC@
IOLIBS=@NETCDF_ALIBS@

models_includes= -I $(top_builddir)/xptor/models \
  -I $(top_builddir)/tglf/src \
  -I $(top_builddir)/neo/src 


#if HAVE_PETSC
#if HAVE_SLEPC
##libEIGEN_la_SOURCES += EIGEN_do.F
#libEIGEN_la_SOURCES += EIGEN_do.F90
#endif
#else
#libEIGEN_la_SOURCES += EIGEN_do.stub.f90
#endif
#
#	
## To get to work with optimization reducing compiler wrappers
## or simply to get .f to work with new automake
AM_LIBTOOLFLAGS = --tag=FC
AM_FFLAGS=$(FC_FIXED_FLAG) $(FC_DBL_FLAG) $(FC_PIC_FLAG) $(models_includes) $(IOINC) $(IONLIBS)
AM_FCFLAGS=$(FC_FREE_FLAG) $(FC_DBL_FLAG) $(FC_PIC_FLAG) $(models_includes) $(IOINC) $(IONLIBS)

CLEANFILES = *.mod *.o *.lo *.lst *.la
