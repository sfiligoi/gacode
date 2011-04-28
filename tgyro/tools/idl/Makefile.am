## ############################################################
##
## File:  Makefile.am
##
## Purpose:  Run automake in top directory to make Makefile.in
##           in all subdirectories
##
## $Id: Makefile.am 50 2010-01-13 05:44:27Z srinath $

## ############################################################


dist_data_DATA = \
	get_tgyro_loc_data.pro \
    	plot_tgyro_loc_summary.pro \
	write_tgyro2gyro_input.pro \
	INPUT.template	
all-local:
	for i in $(dist_noinst_DATA); \
        do \
                cd $(top_builddir)/tgyro/tools/idl && \
                if test ! -f $${i}; \
                  then ln -s $(top_srcdir)/tgyro/tools/idl/$${i} $${i}; \
                fi; \
        done;

MAINTAINERCLEANFILES    = Makefile.in 

   
