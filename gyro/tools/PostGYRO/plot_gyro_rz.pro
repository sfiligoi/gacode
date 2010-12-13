PRO plot_gyro_RZ, data, IT = it, SF=sf, TITLE = title, PLOT_N = plot_n, $
                  PLOT_E = plot_e, PLOT_T = plot_t, $
                  NU_SILENT = nu_silent, _EXTRA = extra
;
; C. Holland, UCSD
; v1.0: 2/27/2007
;
; Given output from GET_GYRO_DATA(), plots fluctuations in (R,Z) plane
; IT specifies timestep, SF specifies interpolation
;
; NOTE: currently only plots phi including zonal flows
;
; v2.0: 3/30/2007: updated to plot density and energy moments; still
; need to remove zonal flows
; v2.1: 5/21/2007: fix for n_theta_plot=1 data
; v3.0: 9/24/07: added plot_t, T = 2/3E - n
; v4.0: 4/8/10: updated to use temp fluctuations from data structure

  DEFAULT, it, data.n_time/2
  DEFAULT, sf, 8

  IF (KEYWORD_SET(plot_t) AND data.exists_t) THEN BEGIN
      field = REFORM(data.mom_t[*,*,data.i_kin,*,it])
      def_title = data.kin_tags[data.i_kin] + ' !4d!XT(R,Z)'
  ENDIF ELSE IF (KEYWORD_SET(plot_E) AND data.exists_e) THEN BEGIN
      field = REFORM(data.mom_e[*,*,data.i_kin,*,it])
      def_title = data.kin_tags[data.i_kin] + ' !4d!XE(R,Z)'
  ENDIF ELSE IF (KEYWORD_SET(plot_n) AND data.exists_n) THEN BEGIN
      field = REFORM(data.mom_n[*,*,data.i_kin,*,it])
      def_title = data.kin_tags[data.i_kin] + ' !4d!Xn(R,Z)'
  ENDIF ELSE BEGIN
      field = REFORM(data.phi[*,*,*,it])
      def_title = '!4u!X(R,Z)'
  ENDELSE

  ;5/21/07: fix here to correct for n_theta_plot = 1 + REFORM above
  IF (data.n_theta_plot EQ 1) THEN BEGIN
      tmpfield = FLTARR(1,data.n_r,data.n_n)
      tmpfield[0,*,*] = field
      field = tmpfield
  ENDIF

  field_rtheta = GYRO_RTHETA_TRANSFORM(field,data.profile_data,sf,$
                                      NU_SILENT=nu_silent)
  DEFAULT, title, def_title + ', it = ' + NUMTOSTRING(it) + ', SF = ' + $
                 NUMTOSTRING(sf)

  GYRO_GENERATE_RZCOORDS, data.profile_data, sf*data.n_theta_plot, R, Z
  GYRO_RZ_COLOR_CONTOUR, field_rtheta, R, Z, TITLE = title, _EXTRA = extra

END ;plot_gyro_RZ

