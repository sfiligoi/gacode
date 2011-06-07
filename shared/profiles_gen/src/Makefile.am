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

noinst_LTLIBRARIES = libgacodprgen_SRC.la
bin_PROGRAMS = prgen

libgacodprgen_SRC_la_SOURCES = \
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
  fluxfit_f_model.f90 
 
prgen_SOURCES = prgen.f90
prgen_LDADD  = libgacodprgen_SRC.la $(top_builddir)/shared/math/libgacodemath.la \
	@FFTW_LIBS@  @LAPACK_LIBS@  @BLAS_LIBS@ $(IOLIBS) 


IOINC=@HDF5_FC_INC@ @NETCDF_INC@
IOLIBS=@NETCDF_ALIBS@ @HDF5_ALIBS@ 

AM_LIBTOOLFLAGS = --tag=FC
AM_FCFLAGS=$(FC_FREE_FLAG) $(FC_DBL_FLAG) $(FC_PIC_FLAG) $(IOINC) -I.
AM_FFLAGS=$(FC_FIXED_FLAG) $(FC_DBL_FLAG) $(FC_PIC_FLAG) $(IOINC) -I.


CLEANFILES = *.mod *.f90 *.lo *.lst *.la


