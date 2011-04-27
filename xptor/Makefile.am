## ############################################################
##
## File:  Makefile.am
##
## Purpose:  Run automake in top directory to make Makefile.in
##           in all subdirectories
##
## $Id: Makefile.am 45 2010-01-12 04:25:07Z srinath $

## ############################################################

SUBDIRS= \
  datread \
  misc \
  inc \
  models

#  idl \
#  nclmisc \
#  ptmisc \
#  solver \
#  source
#  driver \
#libdir = $(prefix)/lib
#lib_LTLIBRARIES = libgacodexptor.la
#
#libgacodexptor_la_SOURCES =
#libgacodexptor_la_LIBADD = src/libgacodexptor_src.la

#if HAVE_MPI
#bin_PROGRAMS     = tglf
#tglf_SOURCES 	= src/tglf.f90
#tglf_LDADD  	= libgacodetglf.la \
#		  @LAPACK_LIBS@ @BLAS_LIBS@
#endif

AM_LIBTOOLFLAGS = --tag=FC
AM_FFLAGS=$(FC_FIXED_FLAG) $(FC_DBL_FLAG)
AM_FCFLAGS=$(FC_FREE_FLAG) $(FC_DBL_FLAG) $(FC_PIC_FLAG)
MAINTAINERCLEANFILES    = Makefile.in
CLEANFILES = *.mod *.o *.lo *.lst *.la
