## ############################################################
##
## File:  Makefile.am
##
## Purpose:  Run automake in top directory to make Makefile.in
##           in all subdirectories
##
## $Id: Makefile.am 176 2010-07-27 15:25:25Z srinath $

## ############################################################

dist_noinst_DATA = \
  inpt.nml \
  xp0in.nml \
  xpin.nml

all-local:
	for i in $(dist_noinst_DATA); \
        do \
                cd $(top_builddir)/xptor/driver && \
                if test ! -f $${i}; \
                  then ln -s $(top_srcdir)/xptor/driver/$${i} $${i}; \
                fi; \
        done;
	

MAINTAINERCLEANFILES    = Makefile.in 

   
