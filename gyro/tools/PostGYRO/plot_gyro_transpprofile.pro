PRO plot_gyro_transpprofile, simdir, ITMIN=itmin, ITMAX=itmax, OPLOT = oplot, $
  TITLE = title, THICK = thick, CHI_E = plot_chi_e, Dn_e = plot_dne, $
  Dn_i = plot_dni, _EXTRA = extra
;
; C. Holland, UCSD
;
; v1.0: 3-5-2007
;
; Given a valid GYRO simdir, plots chi_i vs. r.  Averaged over itmin:itmax
;

  dirpath = GETENV('GYRO_DIR') + '/sim/' + simdir
  
  IF (READ_GYRO_PROFILE_DATA(simdir, profile_data)) THEN BEGIN
      n_r = profile_data.n_r

      time = READ_GYRO_TIMEVECTOR(dirpath)
      n_time = N_ELEMENTS(time)

      DEFAULT, itmin, n_time/2
      DEFAULT, itmax, n_time-1
      DEFAULT, title, simdir + " t = [" + NUMTOSTRING(itmin) + ":" + $
               NUMTOSTRING(itmax) + "]"
      DEFAULT, thick, 1
      NT = itmax-itmin+1
;      array = FLTARR(profile_data.n_kinetic,profile_data.n_field,2,n_r,n_time)
;      IF (READ_GYRO_ARRAY(array, dirpath+'/diff_i.out')) THEN BEGIN
;          chi = REFORM(array[0,0,1,*,*])
;          IF (itmax NE itmin) THEN chi = TOTAL(chi[*,itmin:itmax],2)/NT $
;          ELSE chi = chi[*,itmin]
          IF KEYWORD_SET(plot_dne) THEN BEGIN
              chi = profile_data.diff_ne_sim
              ytitle = 'D!Dne!N/!4v!X!DgB!N'
          ENDIF ELSE IF KEYWORD_SET(plot_dni) THEN BEGIN
              chi = profile_data.diff_ni_sim
              ytitle = 'D!X!Dni!N/!4v!X!DgB!N'
          ENDIF ELSE IF KEYWORD_SET(plot_chi_e) THEN BEGIN
              chi = profile_data.chi_e_sim
              ytitle = '!4v!X!De!N/!4v!X!DgB!N'
          ENDIF ELSE BEGIN
              chi = profile_data.chi_i_sim
              ytitle = '!4v!X!Di!N/!4v!X!DgB!N'
          ENDELSE

          ;if EM sim, add ES and EM components
          IF (profile_data.n_field EQ 2) THEN chi = TOTAL(chi,1)
          IF (itmax EQ itmin) THEN chi = chi[*,itmin] $
          ELSE chi = TOTAL(chi[*,itmin:itmax],2)/NT
          IF KEYWORD_SET(oplot) THEN OPLOT, profile_data.r, chi, $
            THICK = thick, _EXTRA=extra $
          ELSE PLOT, profile_data.r, chi, /XS, XTITLE = 'r/a', $
                     YTITLE = ytitle, CHARSIZE=1.5, $
                     TITLE = title, THICK=thick, XTHICK=thick, YTHICK=thick,$
                     CHARTHICK=thick, _EXTRA = extra
;      ENDIF ELSE BEGIN
;          PRINT, "Couldn't read " + dirpath + "/diff_i.out"
;      ENDELSE
  ENDIF

END ;plot_gyro_transpprofile
