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
gyro.ALCF_BGP \
gyro.FRANKLIN \
gyro.GFORTRAN_PENRYN \
gyro.KESTREL \
gyro.PG_OPT64_FFTW \
gyro.XT5_SLEPC \
gyro.BANACH_PETSC \
gyro.GFORTRAN_TAU \
gyro.KITE \
gyro.PG_OPT64_MUMPS \
gyro.XT_FAST \
gyro.G95_OSX \
gyro.HOPPER2 \
gyro.LOHAN \
gyro.PPPL \
gyro.CARVER \
gyro.GFORTRAN_CORE2 \
gyro.HOPPER2 \
gyro.LOKI \
gyro.REGATTA \
gyro.DELLA \
gyro.GFORTRAN_OSX \
gyro.IFORT_CORE2 \
gyro.UBGL \
gyro.DELPHI \
gyro.GFORTRAN_OSX_BELLI \
gyro.PETREL \
gyro.XT3 \
gyro.DROP \
gyro.GFORTRAN_P4 \
gyro.IFORT_PENRYN \
gyro.PG_OPT64 \
gyro.XT5


all-local:
	for i in $(dist_noinst_DATA) ; \
	do \
		cd $(top_builddir)/shared/exec/ && if test ! -f $${i}; then ln -s $(top_srcdir)/shared/exec/$${i} $${i} ;fi; \
	done


CLEANFILES =  gacodeinput.py

   
