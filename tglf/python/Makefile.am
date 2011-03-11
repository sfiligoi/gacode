## ############################################################
##
## File:  Makefile.am
##
## Purpose:  Run automake in top directory to make Makefile.in
##           in all subdirectories
##
## $Id: Makefile.am 88 2010-04-12 18:43:49Z srinath $

## ############################################################
MAINTAINERCLEANFILES    = Makefile.in 
dist_noinst_DATA= \
   tglf_dict.py
all-local:
	for i in $(dist_noinst_DATA) ; \
	do \
		cd $(top_builddir)/tglf/python/ && if test ! -f $${i}; then ln -s $(top_srcdir)/tglf/python/$${i} $${i} ;fi; \
	done


CLEANFILES = tglf*

