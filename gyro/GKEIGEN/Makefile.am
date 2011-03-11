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

#gyroInc = -I ../src -I ../math @PETSC_INC@ @SLEPC_INC@
#
libdir = $(prefix)/lib/$(FC_LIBSUBDIR)

noinst_LTLIBRARIES = libgacodegkeigen.la
libgacodegkeigen_la_SOURCES = \
  GKEIGEN_globals.f90 \
  GKEIGEN_do.stub.f90 \
  GKEIGEN_matrix_read.f90 \
  GKEIGEN_matrix_write.f90 \
  GKEIGEN_catch_error.f90

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
AM_FFLAGS=$(FC_FIXED_FLAG) $(FC_DBL_FLAG) $(FC_PIC_FLAG) 
AM_FCFLAGS=$(FC_FREE_FLAG) $(FC_DBL_FLAG) $(FC_PIC_FLAG) 

CLEANFILES = *.mod *.o *.lo *.lst *.la
