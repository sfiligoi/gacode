## ##########################################################################
##
## File:		Makefile.am
##
## Author: Mark Fahey
##
## Purpose:  Process this file with automake to produce Makefile
##
## Version:	$Id: Makefile.am 179 2010-09-06 20:18:59Z kruger $
##
## ############################################################################

noinst_LTLIBRARIES = libgacodetgyro_SRC.la
libgacodetgyro_SRC_la_SOURCES = \
  tgyro_globals.f90 \
  tgyro_allocate_globals.f90 \
  tgyro_catch_error.f90 \
  tgyro_comm_setup.f90 \
  tgyro_comm_sync.f90 \
  tgyro_flux.f90 \
  tgyro_flux_vector.f90 \
  tgyro_flux_vector_dense.f90 \
  tgyro_gyro_map.f90 \
  tgyro_init_profiles.f90 \
  tgyro_iteration_variables.f90 \
  tgyro_iteration_parallel.f90 \
  tgyro_iteration_pppl.f90 \
  tgyro_iteration_serial.f90 \
  tgyro_iteration_standard.f90 \
  tgyro_iteration_driver.f90 \
  tgyro_multi_driver.f90 \
  tgyro_neo_map.f90 \
  tgyro_profile_functions.f90 \
  tgyro_qfm_map.f90 \
  tgyro_quasigrad.f90 \
  tgyro_quasineutral.f90 \
  tgyro_read_input.f90 \
  tgyro_readbc_int.f90 \
  tgyro_readbc_real.f90 \
  tgyro_residual.f90 \
  tgyro_restart.f90 \
  tgyro_source.f90 \
  tgyro_stab_driver.f90 \
  tgyro_target_vector.f90 \
  tgyro_tglf_map.f90 \
  tgyro_volume_int.f90 \
  tgyro_write_data.f90 \
  tgyro_write_input.f90 \
  tgyro_write_intermediate.f90

if HAVE_FFTW
THEFFTW = fftw
else
THEFFTW = stub
endif

libdir = $(prefix)/lib
shared_includes= -I $(top_builddir)/shared/math \
	-I $(top_builddir)/shared/GEO \
	-I $(top_builddir)/shared/UMFPACK \
	-I $(top_builddir)/shared/EXPRO  \
	-I $(top_builddir)/shared/vshdf5

gyro_includes = \
	-I $(top_builddir)/gyro/BLEND \
	-I $(top_builddir)/gyro/SSUB \
	-I $(top_builddir)/gyro/TRANSPOSE \
	-I $(top_builddir)/gyro/src 

neo_includes = \
	-I $(top_builddir)/neo/src 

tglf_includes = \
	-I $(top_builddir)/tglf/src 

tgyro_includes = \
  -I $(top_builddir)/tgyro/math \
  -I $(top_builddir)/tgyro/qfm \
  $(shared_includes) \
  $(gyro_includes) \
  $(neo_includes) \
  $(tglf_includes) 


lib_LTLIBRARIES = libgacodestgyro.la
#To get to work with optimization reducing compiler wrappers
IOINC=@HDF5_FC_INC@ @HDF5_INC@
GACODE_ROOT=$(abs_top_builddir)

bin_PROGRAMS     = tgyro


LMATH=@LAPACK_LIBS@ @BLAS_LIBS@ @Z_LIBS@
if USE_FCIOWRAPPERS
if HAVE_FCIOWRAPPERS
IOLIBS=$(FCIOWRAPPERS_ALIBS) @HDF5_ALIBS@ @NETCDF_ALIBS@ @Z_LIBS@
endif
else
IOLIBS=@HDF5_ALIBS@ @NETCDF_ALIBS@ @Z_LIBS@
endif


libgacodestgyro_la_SOURCES =
libgacodestgyro_la_LIBADD = \
    $(top_builddir)/tgyro/math/libgacodetgyromath.la \
    $(top_builddir)/tgyro/qfm/libgacodetgyroqfm.la \
    libgacodetgyro_SRC.la




tgyro_SOURCES = tgyro_main.f90
tgyro_LDADD = libgacodestgyro.la \
              $(top_builddir)/gyro/libgacodesgyro.la \
              $(top_builddir)/neo/src/libgacodesneo.la \
              $(top_builddir)/neo/nclass/libgacodenclass.la \
              $(top_builddir)/tglf/src/libgacodestglf.la \
              $(top_builddir)/shared/libgacodesshared.la \
              @FFTW_LIBS@ @FFTW_LIBS@ $(LMATH) $(IOLIBS)

AM_LIBTOOLFLAGS = --tag=FC
AM_FCFLAGS=$(FC_FREE_FLAG) $(FC_DBL_FLAG) $(FC_PIC_FLAG) $(IOINC)  -I. $(tgyro_includes)
AM_FFLAGS=$(FC_FIXED_FLAG) $(FC_DBL_FLAG) $(FC_PIC_FLAG)  $(IOINC)  -I. $(tgyro_includes)

CLEANFILES = *.mod *.f90 *.lo *.lst *.la tgyro
