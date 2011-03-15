## ############################################################
##
## File:        Makefile.am
##
## Purpose:     Run automake in top directory to make Makefile.in
##              in all subdirectories
##
## Version:     $Id: Makefile.am 3034 2010-04-05 12:42:03Z cary $
##
## Copyright University of Colorado (2000-2005)
##       and Tech-X Corporation (2002-2005)
## See LICENSE file for conditions of use.
##
## ############################################################

# These utilities are important enough that they're required for
# proper building and/or for post installation or input file
# preprocessing and should be installed in the distribution

SCRIPTS_SUFFIXES = .py .sh

# install these scripts
dist_bin_SCRIPTS =  aeqdsk.py geqdsk.py \
	iterDBsemanticMap.py readIterDB.py string_to_number.py

dist_data_DATA =  g118898.03400 a118898.03400

install-exec-local:
	if test "$(abs_top_srcdir)" != "$(abs_top_builddir)"; then \
	  for i in $(dist_bin_SCRIPTS) $(dist_data_DATA); do \
	    if test ! -f $$i; then \
	      ln -s $(srcdir)/$$i .; \
	    fi; \
	  done; \
	fi


## #######
##
## Other targets
##
## #######

include $(top_srcdir)/config/filelists.am
