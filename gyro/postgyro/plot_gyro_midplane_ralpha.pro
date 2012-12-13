PRO plot_gyro_midplane_ralpha, data, IT=it, SF = sf, THICK=thick, $
  PLOT_N = plot_n, PLOT_E = plot_E, PLOT_T = plot_t, $
  PLOT_KINE = plot_kine, PLOT_APAR=plot_apar, PLOT_V=plot_v, $
  FINITE_n = finite_n, LOCAL_NORM=localnorm, $
  FLUXTUBE_AXES=fluxtube_axes, _EXTRA = extra
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
; v5.0: 9/7/11: updated for gacode, added a_||, velocity plotting,
; fluxtube_axis to plot as x/rho_s, y/rho_s rather than r/a, alpha

  DEFAULT, it, data.n_time/2
  DEFAULT, sf, 8

  title = 'midplane '
  IF KEYWORD_SET(finite_n) THEN title += 'finite-n '
  IF (KEYWORD_SET(plot_kine) AND data.exists_n AND data.exists_phi) THEN BEGIN
      field = REFORM(data.mom_n[data.n_theta_plot/2,*,data.n_kinetic-1,*,it]) - $
              REFORM(data.phi[data.n_theta_plot/2,*,*,it])               
      title += '!4d!Xn!De!N/n!D0!N - e!4u!X/T!De!N'
      local_norm = FLTARR(data.n_r) + 1
  ENDIF ELSE IF (KEYWORD_SET(plot_T) AND data.exists_t) THEN BEGIN
      field = REFORM(data.mom_t[data.n_theta_plot/2,*,data.i_kin,*,it])
      title += data.kin_tags[data.i_kin] + ' !4d!XT(r,!4a!X)/T('
      local_norm = data.T_eq[data.i_kin,*]
  ENDIF ELSE IF (KEYWORD_SET(plot_E) AND data.exists_e) THEN BEGIN
      field = REFORM(data.mom_e[data.n_theta_plot/2,*,data.i_kin,*,it])
      title += data.kin_tags[data.i_kin] + ' !4d!XE(r,!4a!X)/p('
      local_norm = data.n_eq[data.i_kin,*]*data.T_eq[data.i_kin,*]
  ENDIF ELSE IF (KEYWORD_SET(plot_n) AND data.exists_n) THEN BEGIN
      field = REFORM(data.mom_n[data.n_theta_plot/2,*,data.i_kin,*,it])
      title += data.kin_tags[data.i_kin] + ' !4d!Xn(r,!4a!X)/n('
      local_norm = data.n_eq[data.i_kin,*]
  ENDIF ELSE IF (KEYWORD_SET(plot_v) AND data.exists_v) THEN BEGIN
      field = REFORM(data.mom_v[data.n_theta_plot/2,*,data.i_kin,*,it])
      title += data.kin_tags[data.i_kin] + ' !4d!Xv(r,!4a!X)/c!Ds!N('
      local_norm = SQRT(data.T_eq[data.n_kinetic-1,*]) ;c_s \prop sqrt(Te)
  ENDIF ELSE IF (KEYWORD_SET(plot_apar) AND data.exists_apar) THEN BEGIN
      field = REFORM(data.Apar[data.n_theta_plot/2,*,*,it])
      title += '(c!Ds!N/c)!4d!XA!D||!N(r,!4a!X)/Te('
      local_norm = SQRT(data.T_eq[data.n_kinetic-1,*])  ;T_e/c_s \prop sqrt(Te)
  ENDIF ELSE IF (data.exists_phi) THEN BEGIN
      field = REFORM(data.phi[data.n_theta_plot/2,*,*,it])
      title += '!4du!X(r,!4a!X)/Te('
      local_norm = data.T_eq[data.n_kinetic-1,*]
  ENDIF ELSE BEGIN
      MESSAGE, 'no viable fluctutaion data selected!', /INFO
      RETURN
  ENDELSE

  IF KEYWORD_SET(finite_n) THEN field[*,0] = 0

  field_ra = GYRO_RALPHA_MIDPLANE_TRANSFORM(field, ALPHA=alpha, $
                                            SCALE_FACTOR=sf)

  IF (KEYWORD_SET(localnorm) AND NOT(KEYWORD_SET(plot_kine))) THEN BEGIN
      datasize = SIZE(field_ra)
      NY = datasize[2]
      IY = FLTARR(NY) + 1
      field_ra /= REFORM(local_norm) # IY
      title += 'r)'
  ENDIF ELSE IF NOT(KEYWORD_SET(plot_kine)) THEN BEGIN
      field_ra /= local_norm[data.n_r/2]
      title += 'r!D0!N)'
  ENDIF
  title += ', t = ' + NUMTOSTRING(it) + ', SF = ' + NUMTOSTRING(sf)

  DEFAULT, thick, 1
  pmax = MAX(ABS(field_ra))
  Ly = 2*!PI*alpha*(data.R0(data.n_r/2)+$
                    data.r(data.n_r/2))/data.n_dn
  Ly = alpha/data.n_dn
  Lx = data.r
  xtitle = 'r!Dmin!N/a'
  ytitle = '!4a!X/2!4p!X'
  r_bnd = [data.r[data.n_bnd], data.r[data.n_r-data.n_bnd-1]]
  IF KEYWORD_SET(fluxtube_axes) THEN BEGIN
      Lx = (data.r - data.r[data.n_r/2])/data.rho_s
      r_bnd = (r_bnd - data.r[data.n_r/2])/data.rho_s
      xtitle = '(r - r!D0!N)/!4q!X!Ds!N'
      Ly = ABS((Ly/MAX(Ly))*2*!PI*data.r[data.n_r/2]/$
               (data.n_dn*data.q[data.n_r/2]*data.rho_s))
      ytitle = '2!4p!X/!4D!X(k!Dy!N!4q!X!Ds!N)'
  ENDIF
print, r_bnd
  COLOR_CONTOUR, field_ra, X_AXIS = Lx, Y_AXIS = Ly, $
                 XTITLE=xtitle, YTITLE=ytitle, $
                 THICK = thick, XTHICK = thick, YTHICK = thick, $
                 CHARTHICK = thick, TITLE = title, $
                 MIN=-pmax,MAX=pmax, _EXTRA = extra
  OPLOT, r_bnd[0]*[1,1],[0,MAX(Ly)], LINESTYLE=2,$
         THICK=thick
  OPLOT,r_bnd[1]*[1,1],[0,Max(Ly)], LINESTYLE=2,$
         THICK=thick
END ;plot_gyro_midplane_ralpha

