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

noinst_LTLIBRARIES = libgacodexptormisc.la
libgacodexptormisc_la_SOURCES = \
  curvss.f \
  dchfev.f \
  dgbtrf.f \
  dgbtrs.f \
  dpchfe.f \
  dpchim.f \
  iarray_max.f \
  inter_cspl.f \
  lspolysmooth.f \
  part.f \
  pchmisc.f \
  pmisc.f \
  r8tomsqz.F \
  rarray_copy.f \
  rarray_scale.f \
  rarray_search.f \
  rarray_sum.f \
  rarray_xdoty.f \
  rarray_zero.f \
  u_erf.f \
  u_lu_backsub.f \
  u_lu_decomp.f \
  u_zeroin.f \
  v_bicubic.f \
  v_spline.f \
  volint.f \
  volintk.f \
  w_bicubic.f \
  w_lin_integ.f \
  w_lin_interp.f \
  w_lin_interp_r8.f \
  w_poly_interp.f \
  w_spline.f \
  write_mod.f

noinst_HEADERS = \
  f77_dcomplx.h


IOINC=@HDF5_FC_INC@ @HDF5_INC@ @NETCDF_INC@
IOLIBS=@HDF5_ALIBS@ @NETCDF_ALIBS@

misc_includes= -I $(top_builddir)/xptor/misc

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
AM_FFLAGS=$(FC_FIXED_FLAG) $(FC_DBL_FLAG) $(FC_PIC_FLAG) $(misc_includes) $(IOINC) $(IONLIBS)
AM_FCFLAGS=$(FC_FREE_FLAG) $(FC_DBL_FLAG) $(FC_PIC_FLAG) $(misc_includes) $(IOINC) $(IONLIBS)

CLEANFILES = *.mod *.o *.lo *.lst *.la
