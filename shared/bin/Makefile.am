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
  gacode_getversion \
  gacode_lines \
  gacode_setup \
  gacode_setup.tcsh \
  gacode_setversion \
  profile_parse.py \
  profiles_gen

all-local:
	for i in $(dist_noinst_DATA) ; \
	do \
		cd $(top_builddir)/shared/bin/ && if test ! -f $${i}; then ln -s $(top_srcdir)/shared/bin/$${i} $${i} ;fi; \
	done


CLEANFILES =  gacodeinput.py

   
