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
  astra2gyro.py \
  astra_map_gyro \
  sampledata.astra

all-local:
	for i in $(dist_noinst_DATA) ; \
	do \
		cd $(top_builddir)/shared/profiles_gen/astra/ && if test ! -f $${i}; then ln -s $(top_srcdir)/shared/profiles_gen/astra/$${i} $${i} ;fi; \
	done


CLEANFILES =  gyro*
#distclean-local:
#	rm -f 

   
