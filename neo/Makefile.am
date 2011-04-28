## ############################################################
##
## File:  Makefile.am
##
## Purpose:  Run automake in top directory to make Makefile.in
##           in all subdirectories
##
## $Id: Makefile.am 81 2010-03-26 17:50:38Z kruger $

## ############################################################

SUBDIRS	= bin nclass src #sim src tools

datadir = $(prefix)/neo
dist_data_DATA = CHANGES

libdir = $(prefix)/lib
lib_LTLIBRARIES = libgacodeneo.la


libgacodeneo_la_SOURCES =
libgacodeneo_la_LIBADD = src/libgacodeneo_src.la

bin_PROGRAMS     = neo
neo_SOURCES = src/neo.f90

neo_LDADD  = libgacodeneo.la \
  $(top_builddir)/neo/nclass/libgacodenclass.la \
  $(top_builddir)/shared/libgacodeshared.la \
 @FFTW_LIBS@  @LAPACK_LIBS@  @BLAS_LIBS@

AM_LIBTOOLFLAGS = --tag=FC
AM_FCFLAGS=$(FC_FREE_FLAG) $(FC_DBL_FLAG) $(FC_PIC_FLAG) $(FFTW_INC) $(shared_includes) -I. 
AM_FFLAGS=$(FC_FIXED_FLAG) $(FC_DBL_FLAG) $(FC_PIC_FLAG) $(FFTW_INC) $(shared_includes) -I.

MAINTAINERCLEANFILES    = Makefile.in 
CLEANFILES = neo *.la *.o

