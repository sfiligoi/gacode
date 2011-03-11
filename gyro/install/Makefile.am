## ###########################################################
##
## File:  Makefile.am
##
## Purpose:  Run automake in top directory to make Makefile.in
##           in all subdirectories
##
## $Id: Makefile.am 176 2010-07-27 15:25:25Z srinath $

## ############################################################

SUBDIRS= 

dist_noinst_DATA= \
	make.ext.ATLAS \
	make.ext.BASSI \
	make.ext.CAESAR \
	make.ext.DAWSON \
	make.ext.DELLA \
	make.ext.DELPHI \
	make.ext.DROP \
	make.ext.FRANKLIN \
	make.ext.G95 \
	make.ext.G95_OSX \
	make.ext.GFORTRAN_CORE2 \
	make.ext.GFORTRAN_OSX \
	make.ext.GFORTRAN_OSX_BELLI \
	make.ext.GFORTRAN_P4 \
	make.ext.GFORTRAN_PENRYN \
	make.ext.GFORTRAN_PM \
	make.ext.GFORTRAN_TAU \
	make.ext.IFORT_P3 \
	make.ext.IFORT_P4 \
	make.ext.IFORT_PM \
	make.ext.JACQUARD \
	make.ext.KESTREL \
	make.ext.KITE \
	make.ext.LF \
	make.ext.PATH_OPT64 \
	make.ext.PETREL \
	make.ext.PG \
	make.ext.PG_FFTW \
	make.ext.PG_OPT64 \
	make.ext.PG_OPT64_FFTW \
	make.ext.PG_OPT64_MUMPS \
	make.ext.PG_P3 \
	make.ext.PHOENIX \
	make.ext.REGATTA \
	make.ext.SEABORG \
	make.ext.UBGL \
	make.ext.XT3 \
	make.ext.XT_CNL \
	make.ext.XT_FAST \
	make.ext.XT_FAST_INLINE \
	make.inc.ATLAS \
	make.inc.BASSI \
	make.inc.CAESAR \
	make.inc.DAWSON \
	make.inc.DELLA \
	make.inc.DELPHI \
	make.inc.DROP \
	make.inc.FRANKLIN \
	make.inc.G95 \
	make.inc.G95_OSX \
	make.inc.GFORTRAN_CORE2 \
	make.inc.GFORTRAN_OSX \
	make.inc.GFORTRAN_OSX_BELLI \
	make.inc.GFORTRAN_P4 \
	make.inc.GFORTRAN_PENRYN \
	make.inc.GFORTRAN_PM \
	make.inc.GFORTRAN_TAU \
	make.inc.IFORT_P3 \
	make.inc.IFORT_P4 \
	make.inc.IFORT_PM \
	make.inc.IFORT_PENRYN\
	make.inc.JACQUARD \
	make.inc.KESTREL \
	make.inc.KITE \
	make.inc.LF \
	make.inc.LOHAN \
	make.inc.PATH_OPT64 \
	make.inc.PETREL \
	make.inc.PG \
	make.inc.PG_FFTW \
	make.inc.PG_OPT64 \
	make.inc.PG_OPT64_FFTW \
	make.inc.PG_OPT64_MUMPS \
	make.inc.PG_P3 \
	make.inc.PHOENIX \
	make.inc.REGATTA \
	make.inc.SEABORG \
	make.inc.UBGL \
	make.inc.XT3 \
	make.inc.XT_CNL \
	make.inc.XT_FAST \
	make.inc.LOKI \
	make.inc.XT5 \
	make.inc.XT5_TAU \
	make.inc.HOPPER \
	make.inc.PC85 \
	make.ext.XT5_TAU \
	make.inc.XT5_HUGEMEM \
	make.inc.XT_FAST_INLINE \
	make.ext.CARVER \
	make.inc.CARVER \
	make.ext.HOPPER \
	make.ext.LOKI \
	make.ext.PC85 \
	make.ext.PPPL \
	make.inc.PPPL \
	make.ext.XT5 \
	make.ext.XT5_HUGEMEM 




all-local:
	for i in $(dist_noinst_DATA) ; \
	do \
		cd $(top_builddir)/install/ && if test ! -f $${i}; then ln -s $(top_srcdir)/install/$${i} $${i} ;fi; \
	done


CLEANFILES = make.inc*

MAINTAINERCLEANFILES    = Makefile.in 

   
