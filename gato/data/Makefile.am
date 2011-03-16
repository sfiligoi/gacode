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
  ingta

all-local:
	for i in $(dist_noinst_DATA) ; \
	do \
		cd $(top_builddir)/gato/data/ && if test ! -f $${i}; then ln -s $(top_srcdir)/gato/data/$${i} $${i} ;fi; \
	done


CLEANFILES = ingta 
#distclean-local:
#	rm -f 

   
