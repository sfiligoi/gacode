## ############################################################
##
## File:  Makefile.am
##
## Purpose:  Run automake in top directory to make Makefile.in
##           in all subdirectories
##
## $Id: Makefile.am 128 2010-06-05 14:47:29Z cary $

## ############################################################

ACLOCAL_AMFLAGS = -I shared/autotool_config

SUBDIRS = shared gato tglf

if HAVE_MPI
SUBDIRS += neo gyro tgyro xptor
endif



MAINTAINERCLEANFILES    = Makefile.in

packagedir = $(prefix)/share
package_DATA = config.summary

install-data-local: $(top_config_summary_file)

distclean-local:
	      rm -f config.summary

uninstall-local:
	      rm -f $(prefix)/$(CONFSHAREDIR)/config.summary
