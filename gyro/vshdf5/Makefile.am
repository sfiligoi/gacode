## ##########################################################################
##
## File:		Makefile.am
##
## Author: Mark Fahey
##
## Purpose:  Process this file with automake to produce Makefile
##
## Version:	$Id: Makefile.am 176 2010-07-27 15:25:25Z srinath $
##
## ############################################################################

noinst_LTLIBRARIES = libgyro_vshdf5.la

gyro_includes= -I $(top_builddir)/gyro/src

# To get to work with optimization reducing compiler wrappers
IOINC=@HDF5_INC@ @FCIOWRAPPERS_INC@
AM_LIBTOOLFLAGS = --tag=FC
AM_FCFLAGS=$(FC_FREE_FLAG) $(FC_DBL_FLAG) $(FC_PIC_FLAG) $(IOINC) $(gyro_includes) $(CPPF90_FLAGS) -I.
AM_FFLAGS=$(FC_FIXED_FLAG) $(FC_DBL_FLAG) $(FC_PIC_FLAG) $(IOINC) $(gyro_includes) $(CPPF90_FLAGS) -I.

if HAVE_FCIOWRAPPERS
libgyro_vshdf5_la_SOURCES = \
	gyro_vshdf5_mod.f90 \
	gyro_write_master_hdf5.f90 
else
libgyro_vshdf5_la_SOURCES = \
	write_hdf5_dummy.f90
endif

CLEANFILES = *.mod *.f90 *.lo *.lst *.la
