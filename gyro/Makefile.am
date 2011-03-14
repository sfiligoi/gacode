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

SUBDIRS	= BLEND TRANSPOSE SSUB src install peek tools bin vugyro vshdf5

datadir = $(prefix)/gyro

dist_data_DATA = CHANGES
#SEK: Exports a bash variable?  Seems odd.
GACODE_ROOT=$(abs_top_builddir)

LMATH=@LAPACK_LIBS@ @BLAS_LIBS@ @Z_LIBS@
IOLIBS=@FCIOWRAPPERS_ALIBS@ @HDF5_ALIBS@ @NETCDF_ALIBS@ @Z_LIBS@

gyro_SOURCES = gyro.f90
gyro_LDADD  = $(top_builddir)/gyro/src/libgacodegyro.la \
	@FFTW_LIBS@  $(LMATH) $(IOLIBS)

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
