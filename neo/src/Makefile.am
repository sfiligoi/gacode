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

libdir = $(prefix)/lib
shared_includes= -I $(top_builddir)/shared/math \
	-I $(top_builddir)/shared/EXPRO \
	-I $(top_builddir)/shared/GEO


noinst_LTLIBRARIES = libgacodeneo_src.la
libgacodeneo_src_la_SOURCES = \
  neo_globals.f90 \
  neo_error.f90 \
  neo_profile_exp.f90 \
  neo_allocate_profile.f90 \
  neo_make_profiles.f90 \
  neo_experimental_profiles.f90 \
  neo_map_experimental_profiles.f90 \
  neo_sparse_solve.f90 \
  neo_compute_fcoll.f90 \
  neo_equilibrium.f90 \
  neo_nclass_dr.f90 \
  neo_energy_grid.f90 \
  neo_g_velocitygrids.f90 \
  neo_theory.f90 \
  neo_check.f90 \
  neo_rotation.f90 \
  neo_transport.f90 \
  neo_interface.f90 \
  neo_init.f90 \
  neo_read_input.f90 \
  neo_run.f90 \
  neo_spitzer.f90 \
  neo_do.f90




 

# To get to work with optimization reducing compiler wrappers
#
AM_LIBTOOLFLAGS = --tag=FC
AM_FCFLAGS=$(FC_FREE_FLAG) $(FC_DBL_FLAG) $(FC_PIC_FLAG) $(FFTW_INC) $(shared_includes) 
AM_FFLAGS=$(FC_FIXED_FLAG) $(FC_DBL_FLAG) $(FC_PIC_FLAG) $(FFTW_INC) $(shared_includes) 


if TX_FORTRAN_MODNAMEISCAP
#include_HEADERS = NEO_INTERFACE.$(TX_FORTRAN_MODEXT)
nodist_include_HEADERS = NEO_INTERFACE.$(TX_FORTRAN_MODEXT)
NEO_INTERFACE.$(TX_FORTRAN_MODEXT): neo_interface.f90
else
#include_HEADERS = neo_interface.$(TX_FORTRAN_MODEXT)
nodist_include_HEADERS = neo_interface.$(TX_FORTRAN_MODEXT)
neo_interface.$(TX_FORTRAN_MODEXT): neo_interface.f90
endif

CLEANFILES = *.mod *.f90 *.lo *.lst *.la


