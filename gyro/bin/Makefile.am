## ############################################################
##
## File:  Makefile.am
##
## Purpose:  Run automake in top directory to make Makefile.in
##           in all subdirectories
##
## $Id: Makefile.am 176 2010-07-27 15:25:25Z srinath $

## ############################################################


MAINTAINERCLEANFILES    = Makefile.in

#sharedir=$(prefix)/share
dist_noinst_DATA= 

dist_bin_SCRIPTS= \
  gyro \
  gyro_clean \
  gyro_parse.py \
  gyro_sim_kill \
  gyro_summary \
  vugyro_dir_init \
  gyro2gkv \
  gyro_kill_do \
  gyro_reg \
  gyro_sim_warn \
  gyro_bat \
  gyro_loc \
  gyro_reg_do \
  simdir_test \
  gyro_genDumpInput.py \
  gyro_genInterface2Global.py

##	gyro_kill_do \

all-local:
	for i in $(dist_bin_SCRIPTS) ; \
	do \
		cd $(top_builddir)/gyro/bin/ && if test ! -f $${i}; then ln -s $(top_srcdir)/gyro/bin/$${i} $${i} ;fi; \
	done
	for i in $(dist_noinst_DATA) ; \
	do \
		cd $(top_builddir)/gyro/bin/ && if test ! -f $${i}; then ln -s $(top_srcdir)/gyro/bin/$${i} $${i} ;fi; \
	done


CLEANFILES = HELP* *gyro* *test*



