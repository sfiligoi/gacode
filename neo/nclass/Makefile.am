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

libdir = $(prefix)/lib/
#
dist_noinst_DATA= \
  pamx_mi.inc \
  pamx_ms.inc \
  pamx_mz.inc

noinst_LTLIBRARIES = libgacodenclass.la
libgacodenclass_la_SOURCES = \
  nclass_mod.f \
  rarray_copy.f \
  rarray_sum.f \
  rarray_zero.f \
  u_erf.f \
  u_lu_backsub.f \
  u_lu_decomp.f \
  write_mod.f
#
all-local:
	for i in $(dist_noinst_DATA) ; \
	  do \
		  cd $(top_builddir)/neo/nclass/ && if test ! -f $${i}; then ln -s $(top_srcdir)/neo/nclass/$${i} $${i} ;fi; \
	  done

# To get to work with optimization reducing compiler wrappers
AM_LIBTOOLFLAGS = --tag=FC
AM_FCFLAGS=$(FC_FREE_FLAG) $(FC_DBL_FLAG) $(FC_PIC_FLAG) $(FFTW_INC)
AM_FFLAGS=$(FC_FIXED_FLAG) $(FC_DBL_FLAG) $(FC_PIC_FLAG) $(FFTW_INC)


CLEANFILES = *.mod *.f90 *.lo *.lst *.la


