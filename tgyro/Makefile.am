## ############################################################
##
## File:  Makefile.am
##
## Purpose:  Run automake in top directory to make Makefile.in
##           in all subdirectories
##
## $Id: Makefile.am 42 2010-01-11 22:35:53Z srinath $

## ############################################################

SUBDIRS	= bin math qfm src #tools

EXTRA_DIST = CHANGES

GACODE_ROOT=$(abs_top_builddir)

bin_PROGRAMS     = tgyro
lib_LTLIBRARIES = libgacodetgyro.la


LMATH=@LAPACK_LIBS@ @BLAS_LIBS@ @Z_LIBS@
if USE_FCIOWRAPPERS
if HAVE_FCIOWRAPPERS
IOLIBS=$(FCIOWRAPPERS_ALIBS) @HDF5_ALIBS@ @NETCDF_ALIBS@ @Z_LIBS@
endif
else
IOLIBS=@HDF5_ALIBS@ @NETCDF_ALIBS@ @Z_LIBS@
endif


libgacodetgyro_la_SOURCES =
libgacodetgyro_la_LIBADD = \
    $(top_builddir)/tgyro/math/libgacodetgyromath.la \
    $(top_builddir)/tgyro/qfm/libgacodetgyroqfm.la \
    $(top_builddir)/tgyro/src/libgacodetgyro_SRC.la

shared_includes= -I $(top_builddir)/shared/math \
	-I $(top_builddir)/shared/GEO \
	-I $(top_builddir)/shared/UMFPACK \
	-I $(top_builddir)/shared/EXPRO  \
	-I $(top_builddir)/shared/vshdf5

gyro_includes = \
	-I $(top_builddir)/gyro/BLEND \
	-I $(top_builddir)/gyro/SSUB \
	-I $(top_builddir)/gyro/TRANSPOSE \
	-I $(top_builddir)/gyro/src

neo_includes = \
	-I $(top_builddir)/neo/src

tglf_includes = \
	-I $(top_builddir)/tglf/src

tgyro_includes = \
  -I $(top_builddir)/tgyro/math \
  -I $(top_builddir)/tgyro/qfm \
  -I $(top_builddir)/tgyro/src \
  $(shared_includes) \
  $(gyro_includes) \
  $(neo_includes) \
  $(tglf_includes)



tgyro_SOURCES = src/tgyro_main.f90
tgyro_LDADD = libgacodetgyro.la \
              $(top_builddir)/gyro/libgacodegyro.la \
              $(top_builddir)/neo/libgacodeneo.la \
              $(top_builddir)/neo/nclass/libgacodenclass.la \
              $(top_builddir)/tglf/libgacodetglf.la \
              $(top_builddir)/shared/libgacodeshared.la \
              @FFTW_LIBS@ @FFTW_LIBS@ $(LMATH) $(IOLIBS)


#To get to work with optimization reducing compiler wrappers
AM_LIBTOOLFLAGS = --tag=FC
AM_FCFLAGS=$(FC_FREE_FLAG) $(FC_DBL_FLAG) $(FC_PIC_FLAG) -I. $(tgyro_includes)
AM_FFLAGS=$(FC_FIXED_FLAG) $(FC_DBL_FLAG) $(FC_PIC_FLAG) -I. $(tgyro_includes)

all-local:
	cd $(top_srcdir)/tgyro/src  && if test ! -h tgyro_main ; then ln -s  $(top_builddir)/tgyro/tgyro tgyro_main ; fi;

MAINTAINERCLEANFILES    = Makefile.in
CLEANFILES = tgyro *.la *.o




