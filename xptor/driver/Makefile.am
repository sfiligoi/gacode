## ############################################################
##
## File:  Makefile.am
##
## Purpose:  Run automake in top directory to make Makefile.in
##           in all subdirectories
##
## $Id: Makefile.am 176 2010-07-27 15:25:25Z srinath $

## ############################################################

dist_noinst_DATA = \
  inpt.nml \
  xp0in.nml \
  xpin.nml

all-local:
	for i in $(dist_noinst_DATA); \
        do \
                cd $(top_builddir)/xptor/driver && \
                if test ! -f $${i}; \
                  then ln -s $(top_srcdir)/xptor/driver/$${i} $${i}; \
                fi; \
        done;
	
libdir = $(prefix)/lib
lib_LTLIBRARIES = libgacodesxptor.la

libgacodesxptor_la_SOURCES =
libgacodesxptor_la_LIBADD = \
  $(top_builddir)/xptor/datread/libgacodexptordatread.la \
  $(top_builddir)/xptor/misc/libgacodexptormisc.la \
  $(top_builddir)/xptor/models/libgacodexptormodels.la \
  $(top_builddir)/xptor/nclmisc/libgacodexptornclmisc.la \
  $(top_builddir)/xptor/ptmisc/libgacodexptorptmisc.la \
  $(top_builddir)/xptor/solver/libgacodexptorsolver.la \
  $(top_builddir)/xptor/source/libgacodexptor_SRC.la

bin_PROGRAMS     = xptor
xptor_SOURCES 	= ptor_driver.f
xptor_LDADD  	= libgacodesxptor.la \
  $(top_builddir)/tglf/libgacodestglf.la \
  $(top_builddir)/neo/libgacodesneo.la \
  $(top_builddir)/shared/libgacodesshared.la \
  @NETCDF_ALIBS@ \
  @LAPACK_LIBS@ @BLAS_LIBS@ 

AM_LIBTOOLFLAGS = --tag=FC
AM_FFLAGS=$(FC_FIXED_FLAG) $(FC_DBL_FLAG)
AM_FCFLAGS=$(FC_FREE_FLAG) $(FC_DBL_FLAG) $(FC_PIC_FLAG)

MAINTAINERCLEANFILES    = Makefile.in 
CLEANFILES = *.mod *.o *.lo *.lst *.la xptor

   
