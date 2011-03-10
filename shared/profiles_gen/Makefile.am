## ############################################################
##
## File:  Makefile.am
##
## Purpose:  Run automake in top directory to make Makefile.in
##           in all subdirectories
##
## $Id: Makefile.am 176 2010-07-27 15:25:25Z srinath $

## ############################################################

SUBDIRS	= BLEND math GEO EXPRO TRANSPOSE SSUB UMFPACK EIGEN src peek tools sim python bin

datadir = $(prefix)/gyro
dist_data_DATA = VERSION

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




