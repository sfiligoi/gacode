## ############################################################
##
## File:  Makefile.am
##
## Purpose:  Run automake in top directory to make Makefile.in
##           in all subdirectories
##
## $Id: Makefile.am 176 2010-07-27 15:25:25Z srinath $

## ############################################################


MAINTAINERCLEANFILES    = Makefile.in 
#
#if HAVE_NETCDF
#bin_PROGRAMS = main
#
#main_SOURCES = \
#	extra.f90
#	read_all.f90
#	main.f90
#
#main_LDADD = $(NETCDF_LIBS)
#endif
## To get to work with optimization reducing compiler wrappers
#AM_LIBTOOLFLAGS = --tag=FC
#AM_FFLAGS = $(FC_FIXED_FLAG) $(INCLUDES)
#AM_FCFLAGS = $(FC_FREE_FLAG) $(INCLUDES)

EXTRA_DIST= \
	extra.f90 \
	read_all.f90 \
	gyro2gkv_exec.f90 \
	gyro2gkv_globals.f90
	main.f90



