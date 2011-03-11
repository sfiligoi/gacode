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
  make.ext.ALCF_BGP \
  make.ext.FRANKLIN  \
  make.ext.HOPPER \
  make.ext.IFORT_PENRYN \
  make.ext.XT5_SLEPC \
  make.ext.DELPHI \
  make.ext.GFORTRAN_OSX_BELLI \
  make.ext.HOPPER2 \
  make.ext.LOKI \
  make.ext.XT_FAST \
  make.ext.DROP \
  make.ext.GFORTRAN_PENRYN \
  make.ext.IFORT_CORE2 \
  make.ext.XT5


all-local:
	for i in $(dist_noinst_DATA) ; \
	do \
		cd $(top_builddir)/gato/install/ && if test ! -f $${i}; then ln -s $(top_srcdir)/gato/install/$${i} $${i} ;fi; \
	done


CLEANFILES = make.ext*
#distclean-local:
#	rm -f 

   
