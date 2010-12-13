;
; C. Holland, UCSD
;
; quick routine to plot linear frequency and growth rates from a GYRO
; simulation.  paramaters should be self-explanitory
;
PRO plot_gyro_linfreq, simdir, ITMIN = itmin, ITMAX = itmax, OPLOT=oplot, $
  VARIANCE = plot_var, THICK=thick, XRANGE = xrange, CHARSIZE=cs, $
  _EXTRA = extra

  err = READ_GYRO_PROFILE_DATA(simdir, profile_data)

  n_n = profile_data.n_n
  ktheta = profile_data.ktheta
  ktheta_label = 'k!D!4h!Nq!X!Ds!N'

  dirpath = GETENV('GYRO_DIR') + '/sim/' + simdir

  time = READ_GYRO_TIMEVECTOR(dirpath)
  n_time = N_ELEMENTS(time)

  DEFAULT, n_n, 1
  array = FLTARR(4,n_n,n_time)
  exists_omega = READ_GYRO_ARRAY(array, dirpath+'/freq.out')

  DEFAULT, itmin, n_time/2
  DEFAULT, itmax, n_time-1
  omega = FLTARR(n_n)
  gamma = FLTARR(n_n)
  omega_err = FLTARR(n_n)
  gamma_err = FLTARR(n_n)
  FOR i_n = 0, n_n-1 DO BEGIN
      omega[i_n] = MEAN(array[0,i_n,itmin:itmax])
      omega_err[i_n] = SDOM(array[0,i_n,itmin:itmax])
      IF KEYWORD_SET(plot_var) THEN $
        omega_err[i_n] = SQRT(VARIANCE(array[0,i_n,itmin:itmax]))
      gamma[i_n] = MEAN(array[1,i_n,itmin:itmax])
      gamma_err[i_n] = SDOM(array[1,i_n,itmin:itmax])
      IF KEYWORD_SET(plot_var) THEN $
        gamma_err[i_n] = SQRT(VARIANCE(array[1,i_n,itmin:itmax]))
  ENDFOR

  DEFAULT, thick, 1
  pmax = MAX(ABS(gamma)) > MAX(ABS(omega))
  IF KEYWORD_SET(oplot) THEN oplot, ktheta, gamma, PSYM=-4, $
    THICk = thick, _EXTRA = extra $
  ELSE PLOT, ktheta, gamma, /XS, XTITLE = ktheta_label, PSYM=-4, $
             YRANGE=[-pmax,pmax], TITLE=simdir, THICK = thick, $
             XTHICK = thick, YTHICK=thick, CHARTHICK = thick, $
             XRANGE = xrange, CHARSIZE=cs, _EXTRA = extra
  ERRPLOT, ktheta, gamma-gamma_err,gamma+gamma_err, THICK = thick, $
           _EXTRA = extra
  OPLOT, ktheta, omega, PSYM=-5, THICK = thick, _EXTRA = extra
  ERRPLOT, ktheta, omega-omega_err, omega+omega_err, THICK = thick,$
           _EXTRA = extra
  OPLOT, ktheta, FLTARR(N_ELEMENTS(ktheta)), LINESTYLE=2, THICK = thick

  IF KEYWORD_SET(oplot) THEN BEGIN
      XYOUTS, 0.1*MAX(ktheta), -0.9*pmax, simdir, CHARTHICK = thick, $
             CHARSIZE=cs, _EXTRA = extra
  ENDIF ELSE BEGIN
      PLOTLABEL, 0.1*MAX(ktheta), 0.9*pmax, 0.15*MAX(ktheta), 0.9*pmax, $
                 ' !4c!X', PSYM=-4, THICK = thick, CHARTHICK=thick,$
                 CHARSIZE=cs, _EXTRA = extra
      PLOTLABEL, 0.1*MAX(ktheta), 0.8*pmax, 0.15*MAX(ktheta), 0.8*pmax, $
                 ' !4x!X', PSYM=-5, THICK = thick, CHARTHICK=thick,$
                 CHARSIZE=cs,_EXTRA = extra
  ENDELSE
END ;plot_gyro_linfreq
