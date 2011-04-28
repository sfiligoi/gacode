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
  comfbl.inc \
  comm30.inc \
  comm31.inc \
  comtr7.inc \
  comtr8.inc \
  comtr9.inc \
  data.m \
  diff \
  efitdat.m \
  gks.m \
  glf.m \
  glf_ave_g.m \
  glf_ave_h.m \
  glf_ave_theta.m \
  glf_dimensions.m \
  glf_eigensolver.m \
  glf_eigenvalues.m \
  glf_get_bd.m \
  glf_get_u.m \
  glf_get_v.m \
  glf_gradB_g.m \
  glf_gradB_h.m \
  glf_hermite.m \
  glf_kpar_g.m \
  glf_kpar_h.m \
  glf_spec-663.m \
  glf_spec.m \
  glf_wd_g.m \
  glf_wd_h.m \
  glf_xgrid.m \
  input.m \
  miller_inputs.m \
  model.m \
  neucg.m \
  pamx_mi.inc \
  pamx_ms.inc \
  pamx_mz.inc \
  pamxkm.inc \
  pamxmp.inc \
  pamxmt.inc \
  pamxnc.inc \
  pamxnx.inc \
  ptor.m \
  s_grid.m \
  share.m \
  sharegk.m \
  tport.m \
  ut.m \
  vtrans.m

all-local:
	for i in $(dist_noinst_DATA); \
        do \
                cd $(top_builddir)/xptor/inc && \
                if test ! -f $${i}; \
                  then ln -s $(top_srcdir)/xptor/inc/$${i} $${i}; \
                fi; \
        done;
	

MAINTAINERCLEANFILES    = Makefile.in 

   
