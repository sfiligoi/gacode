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


libdir = $(prefix)/lib

noinst_LTLIBRARIES = libgacodeblend.la
libgacodeblend_la_SOURCES = \
	BLEND_private.f90 \
	BLEND_init.f90 \
	BLEND_do.f90 \
	BLEND_F.f90 \
	BLEND_Fp.f90 \
	BLEND_f2.f90 \
	BLEND_f2p.f90 \
	BLEND_f3.f90 \
	BLEND_f3p.f90 \
	BLEND_f4.f90 \
	BLEND_f4p.f90 \
	BLEND_cleanup.f90
# To get to work with optimization reducing compiler wrappers
# or simply to get .f to work with new automake
AM_LIBTOOLFLAGS = --tag=FC
AM_FFLAGS=$(FC_FIXED_FLAG) $(FC_DBL_FLAG) $(FC_PIC_FLAG)
AM_FCFLAGS=$(FC_FREE_FLAG) $(FC_DBL_FLAG) $(FC_PIC_FLAG)


CLEANFILES = *.mod *.o *.lo *.lst *.la
