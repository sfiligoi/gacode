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
  README \
  gyroMode2Real.py \
  plotGyro.py \
  readGyroH5.py \
  readGyroH5.txt

all-local:
	for i in $(dist_noinst_DATA) ; \
	do \
		cd $(top_builddir)/gyro/tools/PlotGYROpy && if test ! -f $${i}; then ln -s $(top_srcdir)/gyro/tools/PlotGYROpy/$${i} $${i} ;fi; \
	done


CLEANFILES = *.py *.txt

MAINTAINERCLEANFILES    = Makefile.in 

   
