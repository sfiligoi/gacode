## ############################################################
##
## File:  Makefile.am
##
## Purpose:  Run automake in top directory to make Makefile.in
##           in all subdirectories
##
## $Id: Makefile.am 2 2009-09-01 23:34:41Z srinath $

## ############################################################


MAINTAINERCLEANFILES    = Makefile.in 

dist_bin_SCRIPTS= \
  tgyro \
  tgyro_bat \
  tgyro_clean \
  tgyro_code_lines \
  tgyro_dir \
  tgyro_osg_exec \
  tgyro_parse.py \
  tgyro_plot \
  tgyro_sim_kill \
  tgyro_sim_warn
	  
all-local:
	for i in $(dist_bin_SCRIPTS) ; \
	do \
		cd $(top_builddir)/tgyro/bin/ && if test ! -f $${i}; then ln -s $(top_srcdir)/tgyro/bin/$${i} $${i} ;fi; \
	done

CLEANFILES = tgyro*
