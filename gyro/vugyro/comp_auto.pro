
; do need to start up vugyro widget, then exit, because it runs
;  a number of  COMMON BLOCK INITIALIZATION procedures.
.comp vugyro.pro 
; run vugyro to set up common blocks before compiling routines that
; use the COMMONs:
vugyro,/quick_exit

.comp t_error_event.pro
;  COMMON PRIVATE_SPECTRUM_N is not built by vugyro:
.comp spectrum_n_see.pro
.comp spectrum_n_event.pro

.comp gbflux_see.pro
.comp gbflux_i_see.pro
.comp gbflux_n_see.pro
.comp gbflux_tag.pro
.comp gbflux_event.pro
.comp gbflux_i_event.pro
.comp gbflux_i_plot.pro
.comp gbflux_n_event.pro
.comp get_t_string.pro
.comp t_indices.pro
.comp midplane_fluc_see.pro
.comp midplane_fluc_event.pro
.comp midplane_fluc_plot.pro
.comp osc_see.pro
.comp osc_event.pro
.comp osc_plot.pro
.comp spectrum_np_ave_see.pro
.comp spectrum_np_ave_event.pro
.comp spectrum_np_ave_plot.pro

.comp auto_home_def.pro
.comp main_auto.pro
.comp gbflux_auto.pro
.comp gbflux_i_auto_QG.pro
.comp gbflux_n_auto.pro
.comp midplane_fluc_auto.pro
.comp osc_auto.pro
