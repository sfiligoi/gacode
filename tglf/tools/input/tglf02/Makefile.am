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
  input.tglf \
  input.tglf.gen \
  out.tglf.precision \
  out.tglf.run \
  out.tglf.version 


all-local:
	for i in $(dist_noinst_DATA); \
        do \
                cd $(top_builddir)/tglf/tools/input/tglf02 && \
                if test ! -f $${i}; \
                  then ln -s $(top_srcdir)/tglf/tools/input/tglf02/$${i} $${i}; \
                fi; \
        done;
	

MAINTAINERCLEANFILES    = Makefile.in 

   
