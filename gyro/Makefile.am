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
lib_LTLIBRARIES = libgacodesgyro.la

SUBDIRS	= BLEND TRANSPOSE SSUB tools bin vugyro src GKEIGEN vshdf5 peek 
IOINC =@HDF5_INC@ @NETCDF_INC@

if USE_FCIOWRAPPERS
if HAVE_FCIOWRAPPERS
IOLIBS=$(FCIOWRAPPERS_ALIBS) @HDF5_ALIBS@ @NETCDF_ALIBS@ @Z_LIBS@
IOINC +=$(FCIOWRAPPERS_INC)  
endif
else
IOLIBS= @NETCDF_ALIBS@ @HDF5_ALIBS@ @Z_LIBS@
IOINC+ =-I $(top_builddir)/share/vshdf5  
endif


gyro_includes= -I $(top_builddir)/gyro/src

datadir = $(prefix)/gyro

dist_data_DATA = CHANGES
#SEK: Exports a bash variable?  Seems odd.
GACODE_ROOT=$(abs_top_builddir)

LMATH=@LAPACK_LIBS@ @BLAS_LIBS@ 

libgacodesgyro_la_SOURCES =
libgacodesgyro_la_LIBADD = \
	$(top_builddir)/gyro/src/libgacodegyro_SRC.la \
	$(top_builddir)/gyro/BLEND/libgacodeblend.la \
	$(top_builddir)/gyro/TRANSPOSE/libgacodetranspose.la \
	$(top_builddir)/gyro/SSUB/libgacodessub.la \
	$(top_builddir)/gyro/GKEIGEN/libgacodegkeigen.la \
	$(top_builddir)/gyro/vshdf5/libgyro_vshdf5.la


gyro_SOURCES = src/gyro.f90
gyro_LDADD  = libgacodesgyro.la $(top_builddir)/shared/libgacodesshared.la @FFTW_LIBS@ $(LMATH) $(IOLIBS)


AM_LIBTOOLFLAGS = --tag=FC
AM_FCFLAGS=$(FC_FREE_FLAG) $(FC_DBL_FLAG) $(FC_PIC_FLAG) $(IOINC) $(gyro_includes) $(CPPF90_FLAGS) -I.
AM_FFLAGS=$(FC_FIXED_FLAG) $(FC_DBL_FLAG) $(FC_PIC_FLAG) $(IOINC) $(gyro_includes) $(CPPF90_FLAGS) -I.


MAINTAINERCLEANFILES    = Makefile.in 
CLEANFILES = gyro *.la *.o
