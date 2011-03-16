## ###########################################################
##
## File:  Makefile.am
##
## Purpose:  Run automake in top directory to make Makefile.in
##           in all subdirectories
##
## $Id: Makefile.am 176 2010-07-27 15:25:25Z srinath $

## ############################################################

dist_noinst_DATA= \
make.ext.ALCF_BGP \
make.ext.G95 \
make.ext.HOPPER \
make.ext.LOKI \
make.ext.UBGL \
make.ext.BANACH_PETSC \
make.ext.G95_OSX \
make.ext.HOPPER2 \
make.ext.PATH_OPT64 \
make.ext.XT3 \
make.ext.BASSI \
make.ext.GFORTRAN_CORE2 \
make.ext.IFORT_CORE2 \
make.ext.PETREL \
make.ext.XT5 \
make.ext.CARVER \
make.ext.GFORTRAN_OSX \
make.ext.IFORT_P4 \
make.ext.PG_OPT64 \
make.ext.XT5_SLEPC \
make.ext.DELLA \
make.ext.GFORTRAN_OSX_BELLI \
make.ext.IFORT_PENRYN \
make.ext.PG_OPT64_FFTW \
make.ext.XT_FAST \
make.ext.DELPHI \
make.ext.GFORTRAN_P4 \
make.ext.KESTREL \
make.ext.PG_OPT64_MUMPS \
make.ext.DROP \
make.ext.GFORTRAN_PENRYN \
make.ext.KITE \
make.ext.PPPL \
make.ext.FRANKLIN \
make.ext.GFORTRAN_TAU \
make.ext.LOHAN \
make.ext.REGATTA

all-local:
	for i in $(dist_noinst_DATA) ; \
	do \
		cd $(top_builddir)/gyro/install/ && if test ! -f $${i}; then ln -s $(top_srcdir)/gyro/install/$${i} $${i} ;fi; \
	done


CLEANFILES = make.inc*

MAINTAINERCLEANFILES    = Makefile.in 

   
