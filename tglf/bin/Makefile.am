## ############################################################
##
## File:  Makefile.am
##
## Purpose:  Run automake in top directory to make Makefile.in
##           in all subdirectories
##
## $Id: Makefile.am 88 2010-04-12 18:43:49Z srinath $

## ############################################################

dist_bin_SCRIPTS = \
   HELP \
   tglf \
   tglf_parse.py \
   tglf_sim_warn \
   tglf_version_message 
   

all-local:
	for i in $(dist_bin_SCRIPTS) ; \
	do \
		cd $(top_builddir)/tglf/bin/ && if test ! -f $${i}; then ln -s $(top_srcdir)/tglf/bin/$${i} $${i} ;fi; \
	done


MAINTAINERCLEANFILES    = Makefile.in


