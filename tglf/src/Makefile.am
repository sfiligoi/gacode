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


libdir = $(prefix)/lib
noinst_LTLIBRARIES = libgacodetglf_src.la
libgacodetglf_src_la_SOURCES = \
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

EXTRA_DIST = \
  tglf_inout.f90 \
  tglf_setup_geometry.f90 \
  tglf_TM.f90

lib_LTLIBRARIES = libgacodestglf.la

libgacodestglf_la_SOURCES =
libgacodestglf_la_LIBADD = libgacodetglf_src.la

if HAVE_MPI
bin_PROGRAMS     = tglf
tglf_SOURCES 	= tglf.f90
tglf_LDADD  	= libgacodestglf.la \
		  @LAPACK_LIBS@ @BLAS_LIBS@
endif



if TX_FORTRAN_MODNAMEISCAP
nodist_include_HEADERS= TGLF_INTERFACE.$(TX_FORTRAN_MODEXT)
TGLF_INTERFACE.$(TX_FORTRAN_MODEXT): tglf_interface.f90
else
nodist_include_HEADERS= tglf_interface.$(TX_FORTRAN_MODEXT)
tglf_interface.$(TX_FORTRAN_MODEXT): tglf_interface.f90
endif

# To get to work with optimization reducing compiler wrappers
# or simply for directories with mixed F77 and FC
AM_LIBTOOLFLAGS = --tag=FC
AM_FFLAGS=$(FC_FIXED_FLAG) $(FC_DBL_FLAG)
AM_FCFLAGS=$(FC_FREE_FLAG) $(FC_DBL_FLAG) $(FC_PIC_FLAG)
CLEANFILES = *.mod *.o *.lo *.lst *.la tglf
