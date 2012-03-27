PRO plot_gyro_RZ, data, IT = it, SF=sf, TITLE = title, PLOT_N = plot_n, $
                  PLOT_E = plot_e, PLOT_T = plot_t, PLOT_V = plot_v, $
                  PLOT_Apar=PLOT_apar, NU_SILENT = nu_silent, $
                  LOCAL=local, _EXTRA = extra
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
; v5.0: 8.25.11; updated for gacode compatibility
; v5.1: 9.7.11: added local normalization ability

  DEFAULT, it, data.n_time/2
  DEFAULT, sf, data.theta_mult

  IF (KEYWORD_SET(plot_t) AND data.exists_t) THEN BEGIN
      field = REFORM(data.mom_t[*,*,data.i_kin,*,it])
      def_title = data.kin_tags[data.i_kin] + ' !4d!XT(R,Z)/T!D0!N(r'
      local_norm = data.T_eq[data.i_kin,*]
  ENDIF ELSE IF (KEYWORD_SET(plot_E) AND data.exists_e) THEN BEGIN
      field = REFORM(data.mom_e[*,*,data.i_kin,*,it])
      def_title = data.kin_tags[data.i_kin] + ' !4d!XE(R,Z)/p!D0!N(r'
      local_norm = data.n_eq[data.i_kin,*]*data.T_eq[data.i_kin,*]
  ENDIF ELSE IF (KEYWORD_SET(plot_n) AND data.exists_n) THEN BEGIN
      field = REFORM(data.mom_n[*,*,data.i_kin,*,it])
      def_title = data.kin_tags[data.i_kin] + ' !4d!Xn(R,Z)/n!D0!N(r'
      local_norm = data.n_eq[data.i_kin,*]
  ENDIF ELSE IF (KEYWORD_SET(plot_v) AND data.exists_v) THEN BEGIN
      field = REFORM(data.mom_v[*,*,data.i_kin,*,it])
      def_title = data.kin_tags[data.i_kin] + ' !4d!Xv(R,Z)/c!Ds!N(r'
      local_norm = SQRT(data.T_eq[data.n_kinetic-1,*])
  ENDIF ELSE IF (KEYWORD_SET(plot_Apar) AND data.exists_Apar) THEN BEGIN
      field = REFORM(data.Apar[*,*,*,it])
      def_title = '(c!Ds!N/c)e!4d!XA!D||!N!X(R,Z)/Te(r'
      local_norm = SQRT(data.T_eq[data.n_kinetic-1,*])
  ENDIF ELSE BEGIN
      field = REFORM(data.phi[*,*,*,it])
      def_title = 'e!4du!X(R,Z)/Te(r'
      local_norm = data.T_eq[data.n_kinetic-1,*]
  ENDELSE

  ;5/21/07: fix here to correct for n_theta_plot = 1 + REFORM above
  IF (data.n_theta_plot EQ 1) THEN BEGIN
      tmpfield = FLTARR(1,data.n_r,data.n_n)
      tmpfield[0,*,*] = field
      field = tmpfield
  ENDIF

  field_rtheta_tmp = GYRO_RTHETA_TRANSFORM(field,data,sf,$
                                      NU_SILENT=nu_silent)

  n_r = data.n_r
  ny = sf*data.n_theta_plot
  GYRO_GENERATE_RZCOORDS, data, ny, Rtmp, Ztmp
  field_rtheta = FLTARR(ny+1,n_r)
  R = FLTARR(n_r,ny+1)
  Z = FLTARR(n_r,ny+1)
  field_rtheta[0:ny-1,*] = field_rtheta_tmp
  R[*,0:ny-1] = Rtmp
  Z[*,0:ny-1] = Ztmp
  field_rtheta[ny,*] = field_rtheta_tmp[0,*]
  R[*,ny] = Rtmp[*,0]
  Z[*,ny] = Ztmp[*,0]

  IF KEYWORD_SET(local) THEN BEGIN
      IY = FLTARR(ny+1)+1
      field_rtheta /= IY # REFORM(local_norm)
      def_title += ')'
  ENDIF ELSE BEGIN
      field_rtheta /= local_norm[data.n_r/2]
      def_title += '!D0!N)'
  ENDELSE

  DEFAULT, title, def_title + ', it = ' + NUMTOSTRING(it) + ', SF = ' + $
                 NUMTOSTRING(sf)

  GYRO_RZ_COLOR_CONTOUR, field_rtheta, R, Z, TITLE = title, _EXTRA = extra

END ;plot_gyro_RZ

