## ##########################################################################
##
## File:		Makefile.am
##
## Author: 
##
## Purpose:  Process this file with automake to produce Makefile
##
## Version:	$Id: Makefile.am 59 2010-01-26 00:33:08Z srinath $
##
## ############################################################################


libdir = $(prefix)/lib
noinst_LTLIBRARIES = libgacodetranspose.la
libgacodetranspose_la_SOURCES = \
	fTRANSP_GLOBALS.f90 \
	fTRANSP_INIT.f90 \
	fTRANSP_DO.f90 \
	fTRANSP_CLEANUP.f90 \
	rTRANSP_GLOBALS.f90 \
	rTRANSP_INIT.f90 \
	rTRANSP_DO.f90 \
	rTRANSP_CLEANUP.f90

# To get to work with optimization reducing compiler wrappers
AM_LIBTOOLFLAGS = --tag=FC
AM_FFLAGS=$(FC_FIXED_FLAG) $(FC_DBL_FLAG) $(FC_PIC_FLAG)
AM_FCFLAGS=$(FC_FREE_FLAG) $(FC_DBL_FLAG) $(FC_PIC_FLAG)

CLEANFILES = *.mod *.o *.lo *.lst *.la
