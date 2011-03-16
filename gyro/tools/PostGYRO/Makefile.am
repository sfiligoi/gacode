## ###########################################################
##
## File:  Makefile.am
##
## Purpose:  Run automake in top directory to make Makefile.in
##           in all subdirectories
##
## $Id: Makefile.am 176 2010-07-27 15:25:25Z srinath $

## ############################################################

dist_noinst_DATA= \
  color_contour.pro \
  gyro_rtheta_transform.pro \
  plot_gyro_rz.pro \
  colorbar.pro \
  gyro_rz_color_contour.pro \
  plot_gyro_transpprofile.pro \
  PostGYRO_list.txt \
  default.pro \
  make_bigfields.pro \
  plotlabel.pro \
  analyze_synbes.pro \
  envelope.pro \
  make_synbes_array.pro \
  quick_chi_plot.pro \
  analyze_syncece.pro \
  exp_fit.pro \
  make_syncece_array.pro \
  read_gyro_array.pro \
  analyze_synxphase.pro \
  gauss_fit2.pro \
  make_synrefl_array.pro \
  read_gyro_profile_data.pro \
  calc_coherency1d.pro \
  get_gyro_data.pro \
  numtostring.pro \
  read_gyro_profile_data.pro.new \
  calc_corr1d.pro \
  get_iterdb_data.pro \
  plot_gyro_fluxsurfaces.pro \
  read_gyro_timevector.pro \
  calc_fspect1d.pro \
  gyro_RZ_mpeg.pro \
  plot_gyro_linfreq.pro \
  rms.pro \
  calc_rfspect1d.pro \
  gyro_generate_rzcoords.pro \
  plot_gyro_midplane_ralpha.pro \
  sim_version_test.pro \
  calc_rfspect1d_phase_one.pro \
  gyro_quicksum.pro \
  plot_gyro_midplane_rktspect.pro \
  synbatch.pro \
  calc_xphase1d.pro \
  gyro_ralpha_midplane_transform.pro \
  plot_gyro_midplane_rmsflucamps.pro 


all-local:
	for i in $(dist_noinst_DATA) ; \
	do \
		cd $(top_builddir)/gyro/tools/PostGYRO && if test ! -f $${i}; then ln -s $(top_srcdir)/gyro/tools/PostGYRO/$${i} $${i} ;fi; \
	done


CLEANFILES = make.inc*

MAINTAINERCLEANFILES    = Makefile.in 

   
