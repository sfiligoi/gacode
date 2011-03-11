## ############################################################
##
## File:  Makefile.am
##
## Purpose:  Run automake in top directory to make Makefile.in
##           in all subdirectories
##
## $Id: Makefile.am 187 2010-10-18 18:09:13Z srinath $

## ############################################################


MAINTAINERCLEANFILES    = Makefile.in 
dist_noinst_DATA= \
  make.inc.ALCF_BGP \
  make.inc.G95 \
  make.inc.HOPPER \
  make.inc.LOKI \
  make.inc.UBGL \
  make.inc.BANACH_PETSC \
  make.inc.G95_OSX \
  make.inc.HOPPER2 \
  make.inc.PATH_OPT64 \
  make.inc.XT3 \
  make.inc.BASSI \
  make.inc.GFORTRAN_CORE2 \
  make.inc.IFORT_CORE2 \
  make.inc.PETREL \
  make.inc.XT5 \
  make.inc.CARVER \
  make.inc.GFORTRAN_OSX \
  make.inc.IFORT_P4 \
  make.inc.PG_OPT64 \
  make.inc.XT5_SLEPC \
  make.inc.DELLA \
  make.inc.GFORTRAN_OSX_BELLI \
  make.inc.IFORT_PENRYN \
  make.inc.PG_OPT64_FFTW \
  make.inc.XT_FAST \
  make.inc.DELPHI \
  make.inc.GFORTRAN_P4 \
  make.inc.KESTREL \
  make.inc.PG_OPT64_MUMPS \
  make.inc.DROP \
  make.inc.GFORTRAN_PENRYN \
  make.inc.KITE \
  make.inc.PPPL 

all-local:
	for i in $(dist_noinst_DATA) ; \
	do \
		cd $(top_builddir)/shared/install/ && if test ! -f $${i}; then ln -s $(top_srcdir)/shared/install/$${i} $${i} ;fi; \
	done


CLEANFILES =  gacodeinput.py

   
