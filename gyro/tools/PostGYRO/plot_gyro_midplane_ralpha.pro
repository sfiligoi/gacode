PRO plot_gyro_midplane_ralpha, data, IT=it, SF = sf, THICK=thick, $
  PLOT_N = plot_n, PLOT_E = plot_E, PLOT_T = plot_t, $
  PLOT_KINE = plot_kine, FINITE_n = finite_n, $
  LOCAL_NORM=localnorm, _EXTRA = extra
;
; C. Holland, UCSD
; v1.0 2/27/2007
;
; Plots a snapshot of data field in (r,alpha) at theta = 0, given
; results from GET_GYRO_DATA()
;
; IT is the time index to plot
; SF is an interpolation factor
; THICK sets line thickness (for making postscrpt files)
;
; NOTE: currently only plots phi including zonal flows
;
; v2.0: 3/30/2007 added mom_n, mom_e suppport
; v3.0: 9/27/07: added plot_t, T = 2/3E - n
; v3.1: 10/2/07: plot_kine to plot phi - n_e, finite_n
; v4.0: 4/13/10: update plot_T, local norm

  DEFAULT, it, data.n_time/2
  DEFAULT, sf, 8

  title = 'midplane '
  IF KEYWORD_SET(finite_n) THEN title += 'finite-n '
  IF KEYWORD_SET(plot_kine) THEN BEGIN
      field = REFORM(data.mom_n[data.n_theta_plot/2,*,$
                                data.profile_data.n_kinetic-1,*,it]) - $
              REFORM(data.phi[data.n_theta_plot/2,*,*,it])               
      title += '!4d!Xn!De!N/n!D0!N - e!4u!X/T!De!N'
  ENDIF ELSE IF (KEYWORD_SET(plot_T) AND data.exists_t) THEN BEGIN
      field = REFORM(data.mom_t[data.n_theta_plot/2,*,data.i_kin,*,it])
      title += data.kin_tags[data.i_kin] + ' !4d!XT(r,!4a!X)/T('
      local_norm = data.profile_data.T[data.i_kin,*]
  ENDIF ELSE IF (KEYWORD_SET(plot_E) AND data.exists_e) THEN BEGIN
      field = REFORM(data.mom_e[data.n_theta_plot/2,*,data.i_kin,*,it])
      title += data.kin_tags[data.i_kin] + ' !4d!XE(r,!4a!X)/p('
      local_norm = data.profile_data.n[data.i_kin,*]*data.profile_data.T[data.i_kin,*]
  ENDIF ELSE IF (KEYWORD_SET(plot_n) AND data.exists_n) THEN BEGIN
      field = REFORM(data.mom_n[data.n_theta_plot/2,*,data.i_kin,*,it])
      title += data.kin_tags[data.i_kin] + ' !4d!Xn(r,!4a!X)/n('
      local_norm = data.profile_data.n[data.i_kin,*]
  ENDIF ELSE BEGIN
      field = REFORM(data.phi[data.n_theta_plot/2,*,*,it])
      title += '!4du!X(r,!4a!X)/Te('
      local_norm = data.profile_data.T[data.profile_data.n_kinetic-1,*]
  ENDELSE
  IF KEYWORD_SET(finite_n) THEN field[*,0] = 0

  field_ra = GYRO_RALPHA_MIDPLANE_TRANSFORM(field, ALPHA=alpha, $
                                            SCALE_FACTOR=sf)

  IF KEYWORD_SET(localnorm) THEN BEGIN
      datasize = SIZE(field_ra)
      NY = datasize[2]
      IY = FLTARR(NY) + 1
      field_ra /= REFORM(local_norm) # IY
      title += 'r)'
  ENDIF ELSE BEGIN
      field_ra /= local_norm[data.n_r/2]
      title += 'r!D0!N)'
  ENDELSE
  title += ', t = ' + NUMTOSTRING(it) + ', SF = ' + NUMTOSTRING(sf)

  DEFAULT, thick, 1
  pmax = MAX(ABS(field_ra))
  Ly = 2*!PI*alpha*(data.profile_data.R0(data.n_r/2)+$
                    data.profile_data.r(data.n_r/2))/data.profile_data.n_dn
  Ly = alpha/data.profile_data.n_dn
  COLOR_CONTOUR, field_ra, X_AXIS = data.r, Y_AXIS = Ly, $ ;alpha, $
                 XTITLE = 'r!Dmin!N/a', $; YTITLE = '(2!4p!X/!4D!X!Dn!N)R(r!D0!N)/a', $ ;'!7a!X', $
                 YTITLE = '!4a!X/2!4p!X', $
                 THICK = thick, XTHICK = thick, YTHICK = thick, $
                 CHARTHICK = thick, TITLE = title, $
                 MIN=-pmax,MAX=pmax, _EXTRA = extra
  OPLOT, data.r[data.profile_data.n_bnd]*[1,1],[0,1], LINESTYLE=2,$
         THICK=thick
  OPLOT, data.r[data.profile_data.n_r-data.profile_data.n_bnd]*[1,1],[0,1], LINESTYLE=2,$
         THICK=thick
END ;plot_gyro_midplane_ralpha

