## ############################################################
##
## File:  Makefile.am
##
## Purpose:  Run automake in top directory to make Makefile.in
##           in all subdirectories
##
## $Id: Makefile.am 27 2010-01-11 18:30:13Z srinath $

## ############################################################


MAINTAINERCLEANFILES    = Makefile.in 

dist_noinst_DATA= \
HELP  

dist_bin_SCRIPTS= \
 	neo_parse  \
	neo_sim_warn \
	neo \
	neo_sim_generate \
     	neo_version_message
	

all-local:
	for i in $(dist_bin_SCRIPTS) ; \
	do \
		cd $(top_builddir)/neo/bin/ && if test ! -f $${i}; then ln -s $(top_srcdir)/neo/bin/$${i} $${i} ;fi; \
	done
	for i in $(dist_noinst_DATA) ; \
	do \
		cd $(top_builddir)/neo/bin/ && if test ! -f $${i}; then ln -s $(top_srcdir)/neo/bin/$${i} $${i} ;fi; \
	done


CLEANFILES = HELP* *neo* 



