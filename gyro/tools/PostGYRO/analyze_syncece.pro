PRO analyze_syncece, simdir, NFFT=NFFT, FMIN=fmin, FMAX=fmax, $
                     LOCAL=local,RESULT=result,FILE=file
; C. Holland, UCSD
; v2: 7.27.2010 (v1 is undocumented
;
; cacluates statistics for make_syncece_array output, assumign
; output is saved in array named 'results' in file syncece.sav in
; simdir
;
; by default, assumes working on global simulation.  Set local flag to
; average over radial pairs

  dirpath = GETENV('GYRO_DIR') + '/sim/' + simdir + '/'
  RESTORE, dirpath + 'syncece.sav'
  DEFAULT, NFFT, 32
  N_CECE = results.N_CECE

  df = 1./(NFFT*results.a_over_cs*1e-3/results.interp)
  freq = df*FINDGEN(NFFT/2+1)

  gyro_Te_aspect = FLTARR(2*N_CECE, NFFT/2+1)
  gyro_Te_xspect = FLTARR(N_CECE, NFFT/2+1)
  gyro_Te_xspect_err = FLTARR(N_CECE, NFFT/2+1)
  syn_Te_xspect = FLTARR(N_CECE, NFFT/2+1)
  syn_Te_xspect_err = FLTARR(N_CECE, NFFT/2+1)

  FOR i_tf=0, results.n_tor_frac-1 DO FOR ir = 0, results.N_CECE-1 DO BEGIN
      results.gyro_Te[2*ir,i_tf,*] /= results.Te0[2*ir]
      results.gyro_Te[2*ir+1,i_tf,*] /= results.Te0[2*ir+1]
      results.cece[2*ir,i_tf,*] /= results.Te0[2*ir]
      results.cece[2*ir+1,i_tf,*] /= results.Te0[2*ir+1]

      
      gyro_Te_aspect[2*ir,*] += CALC_RFSPECT1D(results.gyro_Te[2*ir,i_tf,*], $
                                       results.gyro_Te[2*ir,i_tf,*], NFFT=NFFT)
      gyro_Te_aspect[2*ir+1,*] += CALC_RFSPECT1D(results.gyro_Te[2*ir+1,i_tf,*], $
                                       results.gyro_Te[2*ir+1,i_tf,*], NFFT=NFFT)

      gyro_Te_xspect[ir,*] += CALC_RFSPECT1D(results.gyro_Te[2*ir,i_tf,*], $
                                       results.gyro_Te[2*ir+1,i_tf,*], $
                                       NFFT=NFFT, ERR=err)
      gyro_Te_xspect_err[ir,*] += err^2

      syn_Te_xspect[ir,*] += CALC_RFSPECT1D(results.cece[2*ir,i_tf,*], $
                                       results.cece[2*ir+1,i_tf,*], $
                                       NFFT=NFFT, ERR=err)
      syn_Te_xspect_err[ir,*] += err^2
  ENDFOR
  NM = results.N_tor_frac
  gyro_Te_aspect /= NM*df
  gyro_Te_xspect /= NM*df
  gyro_Te_xspect_err = SQRT(gyro_Te_xspect_err/(NM*df))
  syn_Te_xspect /= NM*df
  syn_Te_xspect_err = SQRT(syn_Te_xspect_err/(NM*df))
  apow_rms = 0.

  !P.MULTI = [0,0,2]
  PLOT, freq, gyro_Te_xspect[0,*], /XS, YTITLE = '<|T_e(f)|**2> (1/kHz)',charsize=2
  ERRPLOT, freq, gyro_Te_xspect[0,*]-gyro_Te_xspect_err[0,*], $
           gyro_Te_xspect[0,*]+gyro_Te_xspect_err[0,*]
  OPLOT, freq, gyro_Te_xspect[N_CECE-1,*], LINESTYLE=2
  ERRPLOT, freq, gyro_Te_xspect[N_CECE-1,*]-gyro_Te_xspect_err[N_CECE-1,*], $
           gyro_Te_xspect[N_CECE-1,*]+gyro_Te_xspect_err[N_CECE-1,*]

  OPLOT, freq, syn_Te_xspect[0,*], COLOR=100
  ERRPLOT, freq, syn_Te_xspect[0,*] - syn_Te_xspect_err[0,*], $
           syn_Te_xspect[0,*] + syn_Te_xspect_err[0,*], COLOR=100
  OPLOT, freq, syn_Te_xspect[N_CECE-1,*], LINESTYLE=2, COLOR=100
  ERRPLOT, freq, syn_Te_xspect[N_CECE-1,*]-syn_Te_xspect_err[N_CECE-1,*], $
           syn_Te_xspect[N_CECE-1,*]+syn_Te_xspect_err[N_CECE-1,*], COLOR=100

  DEFAULT, fmin, 0.
  DEFAULT, fmax, MAX(freq)
  mind = MIN(ABS(freq-fmin),idx1)
  mind = MIN(ABS(freq-fmax),idx2)

  gyro_Te_rms = FLTARR(N_CECE)
  PRINT, 'df = ', freq[1]
  PRINT, 'fmin = ', freq[idx1]
  PRINT, 'fmax = ', freq[idx2]

  IF KEYWORD_SET(local) THEN BEGIN
      gyro_Te_aspect = TOTAL(gyro_Te_aspect,1)/(2*N_CECE)
      gyro_Te_xspect = TOTAL(gyro_Te_xspect,1)/N_CECE
      gyro_Te_xspect_err = SQRT(TOTAL(gyro_Te_xspect_err^2,1)/N_CECE)
      syn_Te_xspect = TOTAL(syn_Te_xspect,1)/N_CECE
      syn_Te_xspect_err = SQRT(TOTAL(syn_Te_xspect_err^2,1)/N_CECE)
      
      PLOT, freq, gyro_Te_xspect, xtitle = 'freq (kHz)', charsize=2
      ERRPLOT, freq, gyro_Te_xspect-gyro_Te_xspect_err, $
               gyro_Te_xspect+gyro_Te_xspect_err
      OPLOT, freq, syn_Te_xspect,color=100
      ERRPLOT, freq, syn_Te_xspect-syn_Te_xspect_err,$
               syn_Te_xspect+syn_Te_xspect_err, color=100

      apow_RMS = SQRT(INT_TABULATED(freq, gyro_Te_aspect))
      gyro_Te_rms = SQRT(INT_TABULATED(freq[idx1:idx2],gyro_Te_xspect[idx1:idx2]))
      syn_Te_rms = SQRT(INT_TABULATED(freq[idx1:idx2],syn_Te_xspect[idx1:idx2]))
      gyro_Te_rms_err = 0.5*(INT_TABULATED(freq[idx1:idx2],$
                                             gyro_Te_xspect_err[idx1:idx2]))/$
		    	gyro_Te_rms
      syn_Te_rms_err = 0.5*(INT_TABULATED(freq[idx1:idx2],$
                                              syn_Te_xspect_err[idx1:idx2]))/$
			syn_Te_rms

      PRINT, 'delta S_syn/S_syn: ', INT_TABULATED(freq[idx1:idx2],syn_Te_xspect_err[idx1:idx2])/$
			INT_TABULATED(freq[idx1:idx2],syn_Te_xspect[idx1:idx2])
      
      gyro_te_rms *= 100
      syn_te_rms *= 100
      gyro_te_rms_err *= 100
      syn_te_rms_err *= 100
      apow_RMS *= 100
      PRINT, 'RMS unfiltered T_e (%): ', gyro_Te_rms, $
             '   +/- ', gyro_Te_rms_err
      PRINT, 'RMS unfiltered point T_e (%): ', apow_RMS
      PRINT, 'RMS synthetic T_e (%): ', syn_Te_rms, $
             '   +/- ', syn_Te_rms_err
      
  ENDIF ELSE BEGIN
      gyro_Te_rms = FLTARR(N_CECE)
      syn_Te_rms = FLTARR(N_CECE)
      gyro_Te_rms_err = FLTARR(N_CECE)
      syn_Te_rms_err = FLTARR(N_CECE)

      FOR ir=0,N_CECE-1 DO BEGIN
          gyro_Te_rms[ir]=SQRT(INT_TABULATED(freq[idx1:idx2],gyro_Te_xspect[ir,idx1:idx2]))
          syn_Te_rms[ir]=SQRT(INT_TABULATED(freq[idx1:idx2],syn_Te_xspect[ir,idx1:idx2]))
          gyro_Te_rms_err[ir]=0.5*(INT_TABULATED(freq[idx1:idx2],$
                                                 gyro_Te_xspect_err[ir,idx1:idx2]))/$
				   gyro_Te_rms[ir]
          syn_Te_rms_err[ir] = 0.5*(INT_TABULATED(freq[idx1:idx2],$
                                                  syn_Te_xspect_err[ir,idx1:idx2]))/$
				   syn_Te_rms[ir]
          PRINT, 'r_min/a: ', results.cece_rmid[2*ir]
	  print, 'delta S_syn / S_syn: ', int_Tabulated(freq[idx1:idx2],syn_Te_xspect_err[ir,idx1:idx2])/$
	         int_tabulated(freq[idx1:idx2],syn_Te_xspect[ir,idx1:idx2])
          PRINT, 'RMS unfiltered T_e (%): ', gyro_Te_rms[ir]*100, $
                 '   +/- ', gyro_Te_rms_err[ir]*100
          PRINT, 'RMS synthetic T_e (%): ', syn_Te_rms[ir]*100, $
                 '   +/- ', syn_Te_rms_err[ir]*100
      ENDFOR
      gyro_te_rms *= 100
      syn_te_rms *= 100
      gyro_te_rms_err *= 100
      syn_te_rms_err *= 100

      xloc = results.cece_rmid[2*INDGEN(n_cece)]
      PLOT, xloc, gyro_te_rms, psym=4
      ERRPLOT,xloc, gyro_te_rms-gyro_te_rms_err,$
              gyro_te_rms+gyro_te_rms_err
      
      OPLOT, xloc, syn_te_rms, PSYM=2, color=100
      ERRPLOT, xloc, syn_te_rms-syn_te_rms_err, syn_te_rms+syn_te_rms_err, $
               color=100
     
      oplot, xloc, syn_te_rms*sqrt(1. + 2*syn_te_rms_err/syn_te_rms), color=100,linestyle=2
      oplot, xloc, syn_te_rms*sqrt(1. - 2*syn_te_rms_err/syn_te_rms), color=100,linestyle=2
	
  ENDELSE
  !P.MULTI=0
   result = {$
      freq:freq,$
      apow_RMS:apow_RMS,$
      gyro_te_rms:gyro_te_rms,$
      syn_te_rms:syn_te_rms,$
      gyro_te_rms_err:gyro_te_rms_err,$
      syn_te_rms_err:syn_te_rms_err,$
      gyro_Te_xspect:gyro_Te_xspect,$
      gyro_Te_xspect_err:gyro_Te_xspect_err,$
      syn_Te_xspect:syn_Te_xspect,$
      syn_Te_xspect_err:syn_Te_xspect_err}
      

  IF (KEYWORD_SET(local) AND (N_ELEMENTS(file) NE 0)) THEN BEGIN
      OPENW, 1, file
      PRINTF, 1, '#freq unfilt   unfilt_err    syn    syn_err'
      FOR ii=0,NFFT/2 DO PRINTF, 1, freq[ii], $
        gyro_te_xspect[ii], gyro_te_xspect_err[ii], $
        syn_te_xspect[ii], syn_te_xspect_err[ii]
      CLOSE, 1
  ENDIF
END ;analyze_syncece
