PRO analyze_synrefl, simdir, NFFT=NFFT, FMIN=fmin, FMAX=fmax, $
                     LOCAL=local, IMIN = imin, IMAX=imax
; C. Holland, UCSD
; v2: 7.27.2010 (v1 is undocumented
; v2.1: 4.18.2011
;
; cacluates statistics for make_synrefl_array output, assumign
; output is saved in array named 'results' in file synrefl.sav in
; simdir
;
; by default, assumes working on global simulation.  Set local flag to
; average over radial pairs
;
; use imin/imax to select subrange of CECE pairs to averge over

  dirpath = GETENV('GYRO_DIR') + '/sim/' + simdir + '/'
  RESTORE, dirpath + 'synrefl.sav'
  DEFAULT, NFFT, 32
  N_REFL = results.N_REFL

  df = 1./(NFFT*results.a_over_cs*1e-3/results.interp)
  freq = df*FINDGEN(NFFT/2+1)

  gyro_ne_aspect = FLTARR(2*N_REFL, NFFT/2+1)
  gyro_ne_xspect = FLTARR(N_REFL, NFFT/2+1)
  gyro_ne_xspect_err = FLTARR(N_REFL, NFFT/2+1)
  syn_ne_xspect = FLTARR(N_REFL, NFFT/2+1)
  syn_ne_xspect_err = FLTARR(N_REFL, NFFT/2+1)

  DEFAULT, imin, 0
  DEFAULT, imax, results.N_REFL-1
  FOR i_tf=0, results.n_tor_frac-1 DO FOR ir = imin, imax DO BEGIN      
      gyro_ne_aspect[2*ir,*] += CALC_RFSPECT1D(results.gyro_ne[2*ir,i_tf,*], $
                                       results.gyro_ne[2*ir,i_tf,*], NFFT=NFFT)
      gyro_ne_aspect[2*ir+1,*] += CALC_RFSPECT1D(results.gyro_ne[2*ir+1,i_tf,*], $
                                       results.gyro_ne[2*ir+1,i_tf,*], NFFT=NFFT)

      gyro_ne_xspect[ir,*] += CALC_RFSPECT1D(results.gyro_ne[2*ir,i_tf,*], $
                                       results.gyro_ne[2*ir+1,i_tf,*], $
                                       NFFT=NFFT, ERR=err)
      gyro_ne_xspect_err[ir,*] += err^2

      syn_ne_xspect[ir,*] += CALC_RFSPECT1D(results.refl[2*ir,i_tf,*], $
                                       results.refl[2*ir+1,i_tf,*], $
                                       NFFT=NFFT, ERR=err)
      syn_ne_xspect_err[ir,*] += err^2
  ENDFOR
  NM = results.N_tor_frac
  gyro_ne_aspect /= NM*df
  gyro_ne_xspect /= NM*df
  gyro_ne_xspect_err = SQRT(gyro_ne_xspect_err)/(NM*df)
  syn_ne_xspect /= NM*df
  syn_ne_xspect_err = SQRT(syn_ne_xspect_err)/(NM*df)

  !P.MULTI = [0,0,2]
  PLOT, freq, gyro_ne_xspect[imin,*], /XS, YTITLE = '<|T_e(f)|**2> (1/kHz)',charsize=2
  ERRPLOT, freq, gyro_ne_xspect[imin,*]-gyro_ne_xspect_err[imin,*], $
           gyro_ne_xspect[imin,*]+gyro_ne_xspect_err[imin,*]
  OPLOT, freq, gyro_ne_xspect[imax,*], LINESTYLE=2
  ERRPLOT, freq, gyro_ne_xspect[imax,*]-gyro_ne_xspect_err[imin,*], $
           gyro_ne_xspect[imax,*]+gyro_ne_xspect_err[imax,*]

  OPLOT, freq, syn_ne_xspect[imin,*], COLOR=100
  ERRPLOT, freq, syn_ne_xspect[imin,*] - syn_ne_xspect_err[imin,*], $
           syn_ne_xspect[imin,*] + syn_ne_xspect_err[imin,*], COLOR=100
  OPLOT, freq, syn_ne_xspect[imax,*], LINESTYLE=2, COLOR=100
  ERRPLOT, freq, syn_ne_xspect[imax,*]-syn_ne_xspect_err[imax,*], $
           syn_ne_xspect[imax,*]+syn_ne_xspect_err[imax,*], COLOR=100

  DEFAULT, fmin, 0.
  DEFAULT, fmax, MAX(freq)
  mind = MIN(ABS(freq-fmin),idx1)
  mind = MIN(ABS(freq-fmax),idx2)

  gyro_ne_rms = FLTARR(N_REFL)
  PRINT, 'df = ', freq[1]
  PRINT, 'fmin = ', freq[idx1]
  PRINT, 'fmax = ', freq[idx2]

  IF KEYWORD_SET(local) THEN BEGIN
      N_pair = imax-imin+1
      gyro_ne_aspect = TOTAL(gyro_ne_aspect,1)/(2*N_pair)
      gyro_ne_xspect = TOTAL(gyro_ne_xspect,1)/N_pair
      gyro_ne_xspect_err = SQRT(TOTAL(gyro_ne_xspect_err^2,1)/N_pair)
      syn_ne_xspect = TOTAL(syn_ne_xspect,1)/N_pair
      syn_ne_xspect_err = SQRT(TOTAL(syn_ne_xspect_err^2,1)/N_pair)
      
      PLOT, freq, gyro_ne_xspect, xtitle = 'freq (kHz)'
      ERRPLOT, freq, gyro_ne_xspect-gyro_ne_xspect_err, $
               gyro_ne_xspect+gyro_ne_xspect_err
      OPLOT, freq, syn_ne_xspect,color=100
      ERRPLOT, freq, syn_ne_xspect-syn_ne_xspect_err,$
               syn_ne_xspect+syn_ne_xspect_err, color=100

      apow_RMS = SQRT(INT_TABULATED(freq, gyro_ne_aspect))
      gyro_ne_rms = SQRT(INT_TABULATED(freq[idx1:idx2],gyro_ne_xspect[idx1:idx2]))
      syn_ne_rms = SQRT(INT_TABULATED(freq[idx1:idx2],syn_ne_xspect[idx1:idx2]))
      gyro_ne_rms_err = 0.5*(INT_TABULATED(freq[idx1:idx2],$
                                             gyro_ne_xspect_err[idx1:idx2]))/$
		    	gyro_ne_rms
      syn_ne_rms_err = 0.5*(INT_TABULATED(freq[idx1:idx2],$
                                              syn_ne_xspect_err[idx1:idx2]))/$
			syn_ne_rms

      PRINT, 'delta S_syn/S_syn: ', $
             INT_TABULATED(freq[idx1:idx2],syn_ne_xspect_err[idx1:idx2])/$
             INT_TABULATED(freq[idx1:idx2],syn_ne_xspect[idx1:idx2])

      PRINT, 'RMS unfiltered point n_e (%): ', apow_RMS*100
      PRINT, 'RMS unfiltered n_e (%): ', gyro_ne_rms*100, $
             '   +/- ', gyro_ne_rms_err*100
      PRINT, 'RMS synthetic n_e (%): ', syn_ne_rms*100, $
             '   +/- ', syn_ne_rms_err*100
      
  ENDIF ELSE BEGIN
      gyro_ne_rms = FLTARR(N_REFL)
      syn_ne_rms = FLTARR(N_REFL)
      gyro_ne_rms_err = FLTARR(N_REFL)
      syn_ne_rms_err = FLTARR(N_REFL)
      
      FOR ir=0,N_REFL-1 DO BEGIN
          gyro_ne_rms[ir]=SQRT(INT_TABULATED(freq[idx1:idx2],gyro_ne_xspect[ir,idx1:idx2]))
          syn_ne_rms[ir]=SQRT(INT_TABULATED(freq[idx1:idx2],syn_ne_xspect[ir,idx1:idx2]))
          gyro_ne_rms_err[ir]=0.5*(INT_TABULATED(freq[idx1:idx2],$
                                                 gyro_ne_xspect_err[ir,idx1:idx2]))/$
				   gyro_ne_rms[ir]
          syn_ne_rms_err[ir] = 0.5*(INT_TABULATED(freq[idx1:idx2],$
                                                  syn_ne_xspect_err[ir,idx1:idx2]))/$
				   syn_ne_rms[ir]
          PRINT, 'r_min/a: ', results.refl_rmid[2*ir]
	  print, 'delta S_syn / S_syn: ', $
                 int_tabulated(freq[idx1:idx2],syn_ne_xspect_err[ir,idx1:idx2])/$
                 int_tabulated(freq[idx1:idx2],syn_ne_xspect[ir,idx1:idx2])
          PRINT, 'RMS unfiltered n_e (%): ', gyro_ne_rms[ir]*100, $
                 '   +/- ', gyro_ne_rms_err[ir]*100
          PRINT, 'RMS synthetic n_e (%): ', syn_ne_rms[ir]*100, $
                 '   +/- ', syn_ne_rms_err[ir]*100
      ENDFOR
      gyro_ne_rms *= 100
      syn_ne_rms *= 100
      gyro_ne_rms_err *= 100
      syn_ne_rms_err *= 100

      xloc = results.refl_rmid[2*INDGEN(n_refl)]
      PLOT, xloc, gyro_ne_rms, psym=4
      ERRPLOT,xloc, gyro_ne_rms-gyro_ne_rms_err,$
              gyro_ne_rms+gyro_ne_rms_err
      
      OPLOT, xloc, syn_ne_rms, PSYM=2, color=100
      ERRPLOT, xloc, syn_ne_rms-syn_ne_rms_err, syn_ne_rms+syn_ne_rms_err, $
               color=100
     
      oplot, xloc, syn_ne_rms*sqrt(1. + 2*syn_ne_rms_err/syn_ne_rms), color=100,linestyle=2
      oplot, xloc, syn_ne_rms*sqrt(1. - 2*syn_ne_rms_err/syn_ne_rms), color=100,linestyle=2
	
  ENDELSE
  !P.MULTI=0

END ;analyze_syncece
