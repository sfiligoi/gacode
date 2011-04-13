## ############################################################
##
## File:  Makefile.am
##
## Purpose:  Run automake in top directory to make Makefile.in
##           in all subdirectories
##
## $Id: Makefile.am 176 2010-07-27 15:25:25Z srinath $

## ############################################################

SUBDIRS	= autotool_config envvars math GEO UMFPACK EXPRO profiles_gen vshdf5 python bin exec 

lib_LTLIBRARIES = libgacodeshared.la


libgacodeshared_la_SOURCES =
libgacodeshared_la_LIBADD = \
	$(top_builddir)/shared/math/libgacodemath.la \
	$(top_builddir)/shared/GEO/libgacodegeo.la \
	$(top_builddir)/shared/EXPRO/libgacodeexpro.la \
	$(top_builddir)/shared/UMFPACK/libgacodeumfpack.la \
	$(top_builddir)/shared/vshdf5/libvshdf5.la 


MAINTAINERCLEANFILES  = Makefile.in 
CLEANFILES = *.la *.o




