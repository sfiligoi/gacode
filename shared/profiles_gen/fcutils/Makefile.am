## ############################################################
##
## File:	Makefile.am
##
## Purpose:	Run automake in top directory to make Makefile.in
##		in all subdirectories
##
## Version:	$Id: Makefile.am 449 2007-04-03 18:24:08Z mmiah $
##
## Copyright University of Colorado (2000-2005)
##       and Tech-X Corporation (2002-2005)
## See LICENSE file for conditions of use.
##
## ############################################################

# These utilities are important enough that they're required for
# proper building and/or for post installation or input file
# preprocessing and should be installed in the distribution

# Executable scripts
dist_bin_SCRIPTS =					\
	fcCellMoviePlot.py				\
	fcdata.py						\
	gencoreplot.py					\
	readPsiFit.py					\
	fcplot.py						\
	getIterDbData.py					\
	fcaddpaths.py

# These are imported, but not executable
##dist_data_DATA = \

# Make links to scripts and data if not in place build
$(dist_bin_SCRIPTS): scripts-stamp
all-local: scripts-stamp
scripts-stamp:
	rm -f scripts-temp
	touch scripts-temp
	if test $(ABS_TOP_SRCDIR) != $(abs_top_builddir); then \
	  for i in $(dist_bin_SCRIPTS) $(dist_data_DATA); do \
	    if test -f $(srcdir)/$$i; then \
	      rm -f $$i; \
	      ln -s $(srcdir)/$$i .; \
	    else \
	      echo WARNING: $(srcdir)/$$i does not exist.; \
	    fi; \
	  done; \
	fi
	mv scripts-temp $@

clean-local: clean-scripts
clean-scripts:
	if test $(ABS_TOP_SRCDIR) != $(abs_top_builddir); then \
	  rm -f $(dist_bin_SCRIPTS) $(dist_data_DATA); \
	  rm -f scripts-stamp; \
	fi

## #######
##
## Other targets
##
## #######

SUBDIRS = fcueplot fileconvertors

include $(top_srcdir)/config/filelists.am

