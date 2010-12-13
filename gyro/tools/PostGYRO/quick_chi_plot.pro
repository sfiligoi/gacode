PRO quick_chi_plot, simdir, OPLOT=oplot, THICK=thick, $
                    PRINT_MEAN=print_mean, ITMIN=itmin, ITMAX=itmax, $
                    _EXTRA = extra
;
; C. Holland, UCSD
;
; v1.0: 1-29-2007
; v1.1: 9-11-2007: added thick keyword
; v2.0: 5-5-2008: added print_mean, itmin, itmax keywords
;
; Given a valid GYRO simdir, plots chi_i vs. t.  Useful as a quick
; diagnostic on a particular run
;
; If print_mean flag is set, prints out mean value averaged over 
; [itmin:itmax], default values of [n_time/2:n_time-1]
;
  dirpath = GETENV('GYRO_DIR') + '/sim/' + simdir
  DEFAULT, thick, 1
  IF (READ_GYRO_PROFILE_DATA(simdir, profile_data)) THEN BEGIN
      n_r = profile_data.n_r
      n_theta_plot = profile_data.n_theta_plot
      n_n = profile_data.n_n

      time = READ_GYRO_TIMEVECTOR(dirpath)
      n_time = N_ELEMENTS(time)

      array = FLTARR(profile_data.n_kinetic,profile_data.n_field,2,n_time)
      IF (READ_GYRO_ARRAY(array, dirpath+'/diff.out')) THEN BEGIN
          chi = REFORM(array[0,0,1,*])
          
          IF KEYWORD_SET(oplot) THEN OPLOT, time, chi, $
            THICK = thick, _EXTRA = extra $
          ELSE PLOT, time, chi, /XS, /YS, XTITLE = 't (a/C!Ds!N)', $
                     YTITLE = '!4v!X!Di!N', CHARSIZE=1.5, TITLE = simdir, $
                     THICK=thick,XTHICK=thick,YTHICK=thick,CHARTHICK=thick,$
                     _EXTRA=extra

          IF KEYWORD_SET(print_mean) THEN BEGIN
              DEFAULT, itmin, n_elements(time)/2
              DEFAULT, itmax, n_elements(time)-1
              PRINT, 'Mean chi_i over [' + NUMTOSTRING(time[itmin]) + $
                     ':' + NUMTOSTRING(time[itmax]) +'] = ' + $
                     NUMTOSTRING(MEAN(chi[itmin:itmax]))
              PRINT, 'std. dev of chi_i = ' + $
                     NUMTOSTRING(STDDEV(chi[itmin:itmax]))
          ENDIF

      ENDIF ELSE BEGIN
          PRINT, "Couldn't read diff.out"
      ENDELSE
  ENDIF

END ;quick_chi_plot
