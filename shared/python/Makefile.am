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

instdir=$(prefix)/share
dist_inst_DATA= \
 gacodeinput.py 

all-local:
	for i in $(dist_noinst_DATA) ; \
	do \
		cd $(top_builddir)/shared/python/ && if test ! -f $${i}; then ln -s $(top_srcdir)/shared/python/$${i} $${i} ;fi; \
	done


CLEANFILES =  gacodeinput.py

   
