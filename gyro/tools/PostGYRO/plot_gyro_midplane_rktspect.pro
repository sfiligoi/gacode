PRO plot_gyro_midplane_rktspect, data, ITMIN=itmin, ITMAX=itmax,THICK=thick, $
  PLOT_N = plot_n, PLOT_E = plot_E, PLOT_T = plot_t, $
  FINITE_n = finite_n, LOCAL_NORM=localnorm, $
  PLOT_KSTATS=plot_kstats, KSTATS_FILE=kstats_file, _EXTRA = extra
;
; C. Holland, UCSD
; v1.0 4/13/2010
;

  DEFAULT, itmin, data.n_time/2
  DEFAULT, itmax, data.n_time-1
  title = 'midplane log!D10!N <|'
  IF (KEYWORD_SET(plot_T) AND data.exists_t) THEN BEGIN
      field = REFORM(data.mom_t[data.n_theta_plot/2,*,data.i_kin,*,*])
      title += data.kin_tags[data.i_kin] + ' !4d!XT(r,k!D!4h!X!N)/T('
      local_norm = data.profile_data.T[data.i_kin,*]
  ENDIF ELSE IF (KEYWORD_SET(plot_E) AND data.exists_e) THEN BEGIN
      field = REFORM(data.mom_e[data.n_theta_plot/2,*,data.i_kin,*,*])
      title += data.kin_tags[data.i_kin] + ' !4d!XE(r,k!D!4h!X!N)/p('
      local_norm = data.profile_data.n[data.i_kin,*]*data.profile_data.T[data.i_kin,*]
  ENDIF ELSE IF (KEYWORD_SET(plot_n) AND data.exists_n) THEN BEGIN
      field = REFORM(data.mom_n[data.n_theta_plot/2,*,data.i_kin,*,*])
      title += data.kin_tags[data.i_kin] + ' !4d!Xn(r,k!D!4h!X!N)/n('
      local_norm = data.profile_data.n[data.i_kin,*]
  ENDIF ELSE BEGIN
      field = REFORM(data.phi[data.n_theta_plot/2,*,*,*])
      title += '!4du!X(r,k!D!4h!X!N)/T!De!N('
      local_norm = data.profile_data.T[data.profile_data.n_kinetic-1,*]
  ENDELSE

  Sk = TOTAL(ABS(field[*,*,itmin:itmax])^2,3)/(itmax-itmin+1)
      IY = FLTARR(data.n_n) + 1
  IF KEYWORD_SET(localnorm) THEN BEGIN
      Sk /= REFORM(local_norm) # IY
      title += 'r)'
  ENDIF ELSE BEGIN
      Sk /= local_norm[data.n_r/2]
      title += 'r!D0!N)'
  ENDELSE
  title += '|!U2!N>, t=[' + NUMTOSTRING(data.t[itmin],LEN=5)+':'+$
           NUMTOSTRING(data.t[itmax],LEN=5) + ']'

  DEFAULT, thick, 1
  ktrs =  data.profile_data.rho_s*(data.profile_data.q*$
                                   sqrt(REFORM(data.profile_data.T[data.profile_data.n_kinetic-1,*]))/$
                                   (data.profile_data.r*data.profile_data.B_unit)) # data.n
  ktmax = MAX(ktrs)
  COLOR_CONTOUR, ALOG10(Sk > 1e-16), X_AXIS = data.r#IY, Y_AXIS = ktrs, $ 
                 XTITLE = 'r!Dmin!N/a', YTITLE = 'k!D!4h!Nq!X!Ds!N', $
                 THICK = thick, XTHICK = thick, YTHICK = thick, $
                 CHARTHICK = thick, TITLE = title, $
                 _EXTRA = extra
  OPLOT, data.r[data.profile_data.n_bnd]*[1,1],[0,ktmax], LINESTYLE=2,$
         THICK=thick
  OPLOT, data.r[data.profile_data.n_r-data.profile_data.n_bnd]*[1,1],[0,ktmax], LINESTYLE=2,$
         THICK=thick

  n_r = data.n_r
  avg_k = FLTARR(n_r)
  avg_k2 = FLTARR(n_r)
  rms_dk = FLTARR(n_r)
  FOR i_r = 0, n_r - 1 DO BEGIN
      norm = TOTAL(Sk[i_r,1:*])
      avg_k[i_r] = TOTAL(ktrs[i_r,1:*]*Sk[i_r,1:*])/norm
      avg_k2[i_r] = TOTAL(ktrs[i_r,1:*]^2*Sk[i_r,1:*])/norm
      rms_dk[i_r] = SQRT(avg_k2[i_r] - avg_k[i_r]^2)
  ENDFOR

  IF KEYWORD_SET(plot_kstats) THEN BEGIN
      OPLOT, data.r, avg_k, THICK=thick
      OPLOT, data.r, avg_k + rms_dk, LINESTYLE=4, THICK=thick
      OPLOT, data.r, avg_k - rms_dk, LINESTYLE=4, THICK=thick
  ENDIF 

  IF N_ELEMENTS(kstats_file) NE 0 THEN BEGIN
      OPENW, unit, kstats_file, /GET_LUN
      PRINTF, unit, '# '+data.simdir
      PRINTF, unit, '# t= ' + NUMTOSTRING(data.t[itmin]) + ':' + NUMTOSTRING(data.t[itmax])
      PRINTF, unit, '#r_min/a    <k_t rho_s>    <delta ^2>^(1/2)'
      FOR ii=0,n_r-1 DO PRINTF, unit, data.r[ii],avg_k[ii],rms_dk[ii]
      FREE_LUN, unit 
  ENDIF 
END ;plot_gyro_midplane_rktspect

