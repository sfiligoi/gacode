## ##########################################################################
##
## File:		Makefile.am
##
## Author: Mark Fahey
##
## Purpose:  Process this file with automake to produce Makefile
##
## Version:	$Id: Makefile.am 181 2010-09-06 20:51:58Z kruger $
##
## ############################################################################
#if HAVE_FFTW
#THEFFTW = fftw
#else
#THEFFTW = stub
#endif

libdir = $(prefix)/lib/$(FC_LIBSUBDIR)
shared_includes= -I $(top_builddir)/shared/math \
	-I $(top_builddir)/shared/EXPRO \
	-I $(top_builddir)/shared/GEO
#	-I $(top_builddir)/gyro/BLEND \
#	-I $(top_builddir)/gyro/SSUB \
#	-I $(top_builddir)/gyro/TRANSPOSE \
#	-I $(top_builddir)/gyro/UMFPACK \
#	-I $(top_builddir)/gyro/EIGEN


noinst_LTLIBRARIES = libgacodeneo_src.la
lib_LTLIBRARIES = libgacodeneo.la
libgacodeneo_src_la_SOURCES = \
  neo_globals.f90 \
  neo_profile_exp.f90 \
  neo_allocate_profile.f90 \
  neo_make_profiles.f90 \
  neo_experimental_profiles.f90 \
  neo_map_experimental_profiles.f90 \
  neo_sparse_solve.f90 \
  neo_equilibrium.f90 \
  neo_g_velocitygrids.f90 \
  neo_nclass_dr.f90 \
  neo_energy_grid.f90 \
  neo_theory.f90 \
  neo_check.f90 \
  neo_rotation.f90 \
  neo_transport.f90 \
  neo_interface.f90 \
  neo_init.f90 \
  neo_read_input.f90 \
  neo_run.f90 \
  neo_do.f90

libgacodeneo_la_SOURCES =
libgacodeneo_la_LIBADD = libgacodeneo_src.la

bin_PROGRAMS     = neo
neo_SOURCES = neo.f90
neo_LDADD  = \
  $(top_builddir)/neo/src/libgacodeneo.la \
  $(top_builddir)/neo/nclass/libgacodenclass.la \
  $(top_builddir)/shared/EXPRO/libgacodeexpro.la  \
  $(top_builddir)/shared/math/libgacodemath.la \
  $(top_builddir)/shared/GEO/libgacodegeo.la  \
  $(top_builddir)/shared/UMFPACK/libgacodeumfpack.la  \
  @FFTW_LIBS@  @LAPACK_LIBS@  @BLAS_LIBS@
#neo_INCLUDES = \
#	-I $(top_builddir)/neo/src \
#	-I $(top_builddir)/gyro/GEO 

# To get to work with optimization reducing compiler wrappers
AM_LIBTOOLFLAGS = --tag=FC
AM_FCFLAGS=$(FC_FREE_FLAG) $(FC_DBL_FLAG) $(FC_PIC_FLAG) $(FFTW_INC) $(shared_includes) 
AM_FFLAGS=$(FC_FIXED_FLAG) $(FC_DBL_FLAG) $(FC_PIC_FLAG) $(FFTW_INC) $(shared_includes) 

if TX_FORTRAN_MODNAMEISCAP
include_HEADERS = NEO_INTERFACE.$(TX_FORTRAN_MODEXT)
NEO_INTERFACE.$(TX_FORTRAN_MODEXT): neo_interface.f90
else
include_HEADERS = neo_interface.$(TX_FORTRAN_MODEXT)
neo_interface.$(TX_FORTRAN_MODEXT): neo_interface.f90
endif

CLEANFILES = *.mod *.f90 *.lo *.lst *.la


