## ############################################################
##
## File:  Makefile.am
##
## Purpose:  Run automake in top directory to make Makefile.in
##           in all subdirectories
##
## $Id: Makefile.am 81 2010-03-26 17:50:38Z kruger $

## ############################################################

SUBDIRS	= bin nclass src #sim src tools

datadir = $(prefix)/neo
dist_data_DATA = CHANGES

all-local: link-data-local

link-data-local:
	if test $(abs_top_srcdir) != $(abs_top_builddir); then \
		for i in $(dist_data_DATA); do \
			if test ! -f $$i; then \
				ln -s $(srcdir)/$$i .; \
			fi; \
		done; \
	fi


MAINTAINERCLEANFILES    = Makefile.in 




