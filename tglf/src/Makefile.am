## ##########################################################################
##
## File:		Makefile.am
##
## Author: Mark Fahey
##
## Purpose:  Process this file with automake to produce Makefile
##
## Version:	$Id: Makefile.am 103 2010-04-28 20:50:35Z srinath $
##
## ############################################################################

if HAVE_MPI
bin_PROGRAMS     = tglf
endif

libdir = $(prefix)/lib
lib_LTLIBRARIES = libgacodetglf.la
libgacodetglf_la_SOURCES = \
  tglf_modules.f90 \
  tglf_LS.f90 \
  tglf_geometry.f90 \
  tglf_matrix.f90 \
  tglf_max.f90 \
  tglf_pkg.f90 \
  tglf_interface.f90 \
  tglf_init.f90 \
  tglf_run.f90 \
  tglf_read_input.f90 \
  tglf_dump_input.f90


if HAVE_MPI
tglf_SOURCES 	= tglf.f90
tglf_LDADD  	= $(top_builddir)/tglf/src/libgacodetglf.la \
		  @LAPACK_LIBS@ @BLAS_LIBS@
endif
		 
if TX_FORTRAN_MODNAMEISCAP
include_HEADERS= TGLF_INTERFACE.$(TX_FORTRAN_MODEXT)
TGLF_INTERFACE.$(TX_FORTRAN_MODEXT): tglf_interface.f90
else
include_HEADERS= tglf_interface.$(TX_FORTRAN_MODEXT)
tglf_interface.$(TX_FORTRAN_MODEXT): tglf_interface.f90
endif

EXTRA_DIST = 
# To get to work with optimization reducing compiler wrappers
# or simply for directories with mixed F77 and FC
AM_LIBTOOLFLAGS = --tag=FC
AM_FFLAGS=$(FC_FIXED_FLAG) $(FC_DBL_FLAG)
AM_FCFLAGS=$(FC_FREE_FLAG) $(FC_DBL_FLAG) $(FC_PIC_FLAG)
CLEANFILES = *.mod *.o *.lo *.lst *.la
