## ##########################################################################
##
## File:		Makefile.am
##
## Author: 
##
## Purpose:  Process this file with automake to produce Makefile
##
## Version:	$Id: Makefile.am 176 2010-07-27 15:25:25Z srinath $
##
## ############################################################################


libdir = $(prefix)/lib/$(FC_LIBSUBDIR)
noinst_LTLIBRARIES = libEXPRO.la
libEXPRO_la_SOURCES = \
      	  EXPRO_interface.f90 \
          EXPRO_alloc.f90 \
          EXPRO_alloc_control.f90 \
          EXPRO_palloc.f90 \
          EXPRO_read.f90 \
	  			EXPRO_write_derived.f90 \
					EXPRO_read_driver.f90 \
					EXPRO_qn.f90 \
					EXPRO_write_original.f90 \
          EXPRO_pread.f90

# To get to work with optimization reducing compiler wrappers
AM_LIBTOOLFLAGS = --tag=FC
geoInclude = -I$(top_builddir)/gyro/GEO
AM_FFLAGS=$(FC_FIXED_FLAG) $(FC_DBL_FLAG) $(FC_PIC_FLAG)
AM_FCFLAGS=$(FC_FREE_FLAG) $(FC_DBL_FLAG) $(FC_PIC_FLAG) $(geoInclude)
CLEANFILES = *.mod *.o *.lo *.lst *.la
