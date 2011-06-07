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
SUBDIRS = data 

libdir = $(prefix)/lib

AM_LIBTOOLFLAGS = --tag=FC
AM_FCFLAGS=$(FC_FREE_FLAG) $(FC_DBL_FLAG) $(FC_PIC_FLAG) -I.
AM_FFLAGS=$(FC_FIXED_FLAG) $(FC_DBL_FLAG) $(FC_PIC_FLAG) -I.




bin_PROGRAMS     = gato
gato_SOURCES = \
  gato_plotg_library.stub.f \
  gato.f 

CLEANFILES = *.mod *.f90 *.lo *.lst *.la *.f


