## ##########################################################################
##
## File:		Makefile.am
##
## Author: Mark Fahey
##
## Purpose:  Process this file with automake to produce Makefile
##
## Version:	$Id: Makefile.am 180 2010-09-06 20:50:31Z kruger $
##
## ############################################################################

libdir = $(prefix)/lib/$(FC_LIBSUBDIR)
bin_PROGRAMS     = prgen

IOINC=@HDF5_FC_INC@ @FCIOWRAPPERS_INC@ @NETCDF_INC@
IOLIBS=@FCIOWRAPPERS_ALIBS@ @HDF5_ALIBS@ @NETCDF_ALIBS@

AM_LIBTOOLFLAGS = --tag=FC
AM_FCFLAGS=$(FC_FREE_FLAG) $(FC_DBL_FLAG) $(FC_PIC_FLAG) $(IOINC) $(gyro_includes) $(CPPF90_FLAGS) -I.
AM_FFLAGS=$(FC_FIXED_FLAG) $(FC_DBL_FLAG) $(FC_PIC_FLAG) $(IOINC) $(gyro_includes) $(CPPF90_FLAGS) -I.




prgen_SOURCES = \
  prgen_read_globals.f90 \
  prgen_fluxfit_globals.f90 \
  allocate_internals.f90 \
  allocate_iterdb_vars.f90 \
  allocate_plasmastate_vars.f90 \
  allocate_peqdsk_vars.f90 \
  prgen_read_iterdb.f90 \
  prgen_read_iterdb_nc.f90 \
  prgen_map_iterdb.f90 \
  prgen_map_plasmastate.f90 \
  prgen_map_peqdsk.f90 \
  prgen_read_cer.f90 \
  prgen_read_gato.f90 \
  prgen_read_plasmastate.f90 \
  prgen_read_peqdsk.f90 \
  prgen_write.f90 \
  fluxfit_driver.f90 \
  fluxfit_moments.f90 \
  fluxfit_intersect.f90 \
  fluxfit_minmax.f90 \
  fluxfit_fourier.f90 \
  fluxfit_error.f90 \
  fluxfit_f_model.f90 \
  prgen.f90

prgen_LDADD  = $(top_builddir)/shared/math/libgacodesmath.la \
	@FFTW_LIBS@  @LAPACK_LIBS@  @BLAS_LIBS@ $(IOLIBS) \
	$(top_builddir)/gyro/EXPRO/libEXPRO.la

BigScience_INCLUDES = \
	-I $(top_builddir)/gyro/src \
	-I $(top_builddir)/gyro/BLEND \
	-I $(top_builddir)/gyro/math \
	-I $(top_builddir)/gyro/TRANSPOSE \
	-I $(top_builddir)/gyro/SSUB \
	-I $(top_builddir)/gyro/GEO \
	-I $(top_builddir)/gyro/UMFPACK \
	-I $(top_builddir)/gyro/EXPRO

if TX_FORTRAN_MODNAMEISCAP
include_HEADERS = GYRO_GLOBALS.$(TX_FORTRAN_MODEXT) GYRO_INTERFACE.$(TX_FORTRAN_MODEXT)
GYRO_GLOBALS.$(TX_FORTRAN_MODEXT):   gyro_globals.f90
GYRO_INTERFACE.$(TX_FORTRAN_MODEXT): gyro_interface.f90
else
include_HEADERS = gyro_globals.$(TX_FORTRAN_MODEXT) gyro_interface.$(TX_FORTRAN_MODEXT)
gyro_globals.$(TX_FORTRAN_MODEXT):   gyro_globals.f90
gyro_interface.$(TX_FORTRAN_MODEXT): gyro_interface.f90
endif


CLEANFILES = *.mod *.f90 *.lo *.lst *.la

write_hdf5.o: hdf5_mod.o
write_hdf5.lo: hdf5_mod.lo
gyro_fieldeigen.o: gyro_fieldeigen_private.o
gyro_fieldeigen.lo: gyro_fieldeigen_private.lo
