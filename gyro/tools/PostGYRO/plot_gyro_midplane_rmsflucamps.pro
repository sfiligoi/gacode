PRO plot_gyro_midplane_rmsflucamps, data, ITMIN=itmin, ITMAX=itmax, $
  PLOT_N=plot_n, PLOT_E=plot_e, PLOT_T=plot_T, OPLOT=overplot, $
  ERRPLOT=errplot, LOCALNORM=localnorm, FINITE_N = finite_n, THICK=thick, $
  RAVG_RANGE = ravg_range, FILENAME=filename, IN1=in1, IN2=in2, _EXTRA=extra
;
; C. Holland, UCSD
; v1.0: 7-16-2007
;
; Plots midplane RMS values of n > 0 fluctuations.  Use localnorm to
; plot f_tilde(r)/f0(r) instead of f_tilde(r)/f0(r0)
;
; v2.0: 9-24-07: revised for errors, added plot_t, T = 2/3E - n
;                amp_err def should be checked.  Basic idea:
;                y = <|f|^2> = y_avg +/- eps
;                z = sqrt(y) ~ sqrt(y_avg) +/- eps/2*sqrt(y_avg)
;                so rmsamp_err = eps/2*sqrt(rms_amp); use SDOM for
;                calculting eps.
;
;                Use def: f(theta=alpha=0) = Re{f_n=0 + 2*sum on n
;                f_n>0}
;                Set ravg_range = [rmin,rmax] to define an averaging range
;
; v3.0: 3-16-2010: updated to use correct temperature fluctuations
;
  DEFAULT, itmin, data.n_time/2
  DEFAULT, itmax, data.n_time-1
  n_time = itmax-itmin+1
  DEFAULT, thick, 1
  DEFAULT, in1, 1
  DEFAULT, in2, data.n_n-1

  IF (KEYWORD_SET(plot_T) AND data.exists_t) THEN BEGIN
      field = REFORM(data.mom_t[data.n_theta_plot/2,*,data.i_kin, $
                                    *,itmin:itmax])
      def_title = data.kin_tags[data.i_kin] + ' !4d!XT/T'
      IF KEYWORD_SET(localnorm) THEN BEGIN
          lnorm =data.profile_data.T[data.i_kin,*]
          def_title += '(r)'
      ENDIF ELSE BEGIN
          lnorm = data.profile_data.T[data.i_kin,data.n_r/2]
          def_title += '(r!D0!N)'
      ENDELSE

  ENDIF ELSE IF (KEYWORD_SET(plot_E) AND data.exists_e) THEN BEGIN
      field = REFORM(data.mom_e[data.n_theta_plot/2,*,data.i_kin, $
                                    *,itmin:itmax])
      def_title = data.kin_tags[data.i_kin] + ' (3/2)!4d!Xp/p'
      IF KEYWORD_SET(localnorm) THEN BEGIN
          lnorm = data.profile_data.n[data.i_kin,*]*$
                  data.profile_data.T[data.i_kin,*]
          def_title += '(r)'
      ENDIF ELSE BEGIN
          lnorm = data.profile_data.n[data.i_kin,data.n_r/2]*$
                  data.profile_data.T[data.i_kin,data.n_r/2]
          def_title += '(r!D0!N)'
      ENDELSE
  ENDIF ELSE IF (KEYWORD_SET(plot_n) AND data.exists_n) THEN BEGIN
      field = REFORM(data.mom_n[data.n_theta_plot/2,*,data.i_kin, $
                                    *,itmin:itmax])
      def_title = data.kin_tags[data.i_kin] + ' !4d!Xn/n'
      IF KEYWORD_SET(localnorm) THEN BEGIN
          lnorm = data.profile_data.n[data.i_kin,*]
          def_title += '(r)'
      ENDIF ELSE BEGIN
          lnorm = data.profile_data.n[data.i_kin,data.n_r/2]
          def_title += '(r!D0!N)'
      ENDELSE
  ENDIF ELSE BEGIN
      field = REFORM(data.phi[data.n_theta_plot/2,*,*,itmin:itmax])
      def_title = 'e!4u!X/T!De!N'
      IF KEYWORD_SET(localnorm) THEN BEGIN
          lnorm = data.profile_data.T[data.profile_data.n_spec-1,*]
          def_title += '(r)'
      ENDIF ELSE BEGIN
          lnorm = data.profile_data.T[data.profile_data.n_spec-1,data.n_r/2]
	  def_title += '(r!D0!N)'
      ENDELSE	
  ENDELSE

  pwr = 2*TOTAL(ABS(field[*,in1:in2,*])^2,2)
  IF NOT(KEYWORD_SET(finite_n)) THEN pwr += $
	REFORM(ABS(field[*,0,*])^2)
  amp = SQRT(TOTAL(pwr,2)/(itmax-itmin+1))
  amp_err = 0

  amp /= lnorm
  amp_err /= lnorm
 
  ;change to percent
  amp *= 100
  amp_err *=100

  ;make title  
  IF NOT(KEYWORD_SET(finite_n)) THEN DEFAULT, title, 'RMS ' + def_title $
  ELSE DEFAULT, title, 'RMS finite-n ' + def_title + ' (%)'

  IF KEYWORD_SET(overplot) THEN BEGIN
      OPLOT, data.r, amp, THICK=thick, _EXTRA=extra
  ENDIF ELSE BEGIN
      PLOT, data.r, amp, /XS, XTITLE ='r/a', TITLE=title, $
            THICK=thick,XTHICK=thick,YTHICK=thick,CHARTHICK=thick, $
            YRANGE = [0, MAX(amp+amp_err)], _EXTRA=extra
  ENDELSE
  IF KEYWORD_SET(errplot) THEN BEGIN
      OPLOT, data.r, amp-amp_err, THICK=thick, LINESTYLE=2, _EXTRA=extra
      OPLOT, data.r, amp+amp_err, THICK=thick, LINESTYLE=2, _EXTRA=extra
  ENDIF

  n_bnd = data.profile_data.n_bnd
  IF N_ELEMENTS(ravg_range) EQ 2 THEN ar = ravg_range $
  ELSE ar = [data.r[n_bnd], data.r[data.n_r-1-n_bnd]]
  idx = WHERE((data.r GE ar[0]) AND (data.r LE ar[1]))
  amp_mean = MEAN(amp[idx])
  amp_mean_err = SDOM(amp[idx])
  XYOUTS, data.r[data.n_r/10], 0.2*MAX(amp), "r-avg mean: " + $
          NUMTOSTRING(amp_mean) + " !M+ " + NUMTOSTRING(amp_mean_err), $
          CHARTHICK=thick, _EXTRA=extra
  OPLOT, ar, [1,1]*amp_mean, THICK=thick, _EXTRA=extra
  OPLOT, ar, [1,1]*amp_mean+amp_mean_err, THICK=thick, LINESTYLE=2, $
         _EXTRA=extra
  OPLOT, ar, [1,1]*amp_mean-amp_mean_err, THICK=thick, LINESTYLE=2, $
         _EXTRA=extra

  PRINT, 'mean RMS value = ', amp_mean, ' +/- ', amp_mean_err

  IF N_ELEMENTS(filename) NE 0 THEN BEGIN
	OPENW, 1, GETENV('GYRO_DIR') + '/sim/' + data.simdir + '/' + filename
	PRINTF, 1, data.n_r
	FOR i_r = 0, data.n_r-1 DO PRINTF, 1, data.r[i_r], amp[i_r]
	CLOSE, 1
  ENDIF

END ;plot_gyro_midplane_rmsflucamps
