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
  gyro_plotH5.py \
  plot-gyro3D.py
 
all-local:
	for i in $(dist_noinst_DATA) ; \
	do \
		cd $(top_builddir)/gyro/tools/PlotGYROpy&& if test ! -f $${i}; then ln -s $(top_srcdir)/gyro/tools/PlotGYROpy/$${i} $${i} ;fi; \
	done


CLEANFILES = make.inc*

MAINTAINERCLEANFILES    = Makefile.in 

   
