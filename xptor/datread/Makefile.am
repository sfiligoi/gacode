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

noinst_LTLIBRARIES = libgacodexptordatread.la
libgacodexptordatread_la_SOURCES = \
  average7_1d.f \
  cdfread.f \
  cub_spline.f90 \
  datavg.f \
  datexch.f \
  datmap.f \
  datnclass.f \
  datneo.f \
  gridsetup.f \
  neutxp.f \
  read_er.f \
  readiterdb.f \
  readtestprofiles.f \
  readufiles.f \
  tenspline.f \
  u0read.f \
  u1read.f \
  u1tread.f \
  u2read.f \
  u2tread.f \
  u3read.f

IOINC= @NETCDF_INC@
IOLIBS=@NETCDF_ALIBS@

datread_includes= -I $(top_builddir)/xptor/datread

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
AM_FFLAGS=$(FC_FIXED_FLAG) $(FC_DBL_FLAG) $(FC_PIC_FLAG) $(datread_includes) $(IOINC) $(IONLIBS)
AM_FCFLAGS=$(FC_FREE_FLAG) $(FC_DBL_FLAG) $(FC_PIC_FLAG) $(datread_includes) $(IOINC) $(IONLIBS)

CLEANFILES = *.mod *.o *.lo *.lst *.la
