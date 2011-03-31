## ############################################################
##
## File:  Makefile.am
##
## Purpose:  Run automake in top directory to make Makefile.in
##           in all subdirectories
##
## $Id: Makefile.am 176 2010-07-27 15:25:25Z srinath $

## ############################################################

bin_PROGRAMS=gyro
lib_LTLIBRARIES = libgacodegyro.la

SUBDIRS	= BLEND TRANSPOSE SSUB src tools bin vugyro GKEIGEN vshdf5 peek

gyro_includes= -I $(top_builddir)/gyro/src
AM_LIBTOOLFLAGS = --tag=FC
AM_FCFLAGS=$(FC_FREE_FLAG) $(FC_DBL_FLAG) $(FC_PIC_FLAG) $(IOINC) $(gyro_includes) $(CPPF90_FLAGS) -I.
AM_FFLAGS=$(FC_FIXED_FLAG) $(FC_DBL_FLAG) $(FC_PIC_FLAG) $(IOINC) $(gyro_includes) $(CPPF90_FLAGS) -I.

datadir = $(prefix)/gyro

dist_data_DATA = CHANGES
#SEK: Exports a bash variable?  Seems odd.
GACODE_ROOT=$(abs_top_builddir)

LMATH=@LAPACK_LIBS@ @BLAS_LIBS@ @Z_LIBS@
IOLIBS=@HDF5_ALIBS@ @NETCDF_ALIBS@ @Z_LIBS@

libgacodegyro_la_SOURCES =
libgacodegyro_la_LIBADD = \
	$(top_builddir)/gyro/src/libgacodegyro_SRC.la \
	$(top_builddir)/gyro/BLEND/libgacodeblend.la \
	$(top_builddir)/shared/math/libgacodemath.la \
	$(top_builddir)/gyro/TRANSPOSE/libgacodetranspose.la \
	$(top_builddir)/gyro/SSUB/libgacodessub.la \
	$(top_builddir)/shared/GEO/libgacodegeo.la \
	$(top_builddir)/shared/EXPRO/libgacodeexpro.la \
	$(top_builddir)/shared/UMFPACK/libgacodeumfpack.la \
	$(top_builddir)/gyro/GKEIGEN/libgacodegkeigen.la \
	$(top_builddir)/shared/vshdf5/libvshdf5.la \
	$(top_builddir)/gyro/vshdf5/libgyro_vshdf5.la


gyro_SOURCES = src/gyro.f90
gyro_LDADD  = libgacodegyro.la @FFTW_LIBS@  $(LMATH) $(IOLIBS)

all-local: link-data-local

link-data-local:
	if test $(abs_top_srcdir) != $(abs_top_builddir); then \
		for i in $(dist_data_DATA); do \
			if test ! -f $$i; then \
				ln -s $(srcdir)/$$i .; \
			fi; \
		done; \
	fi

MAINTAINERCLEANFILES    = Makefile.in 
CLEANFILES = gyro *.la *.o
