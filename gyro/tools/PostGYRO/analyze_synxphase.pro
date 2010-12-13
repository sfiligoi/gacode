PRO analyze_synxphase, simdir, NFFT=NFFT, FMIN=fmin, FMAX=fmax, result=result,$
                       FILE=file

; C. Holland, UCSD
; v1.0: 9/17/2010: based of 2009 APS routine by AEW, uses output from
; make_syncece_arrays.pro, make_synrefl_arrays.pro
;

  dirpath = GETENV('GYRO_DIR') + '/sim/' + simdir + '/'
  RESTORE, dirpath + 'syncece.sav'
  cece = results
  RESTORE, dirpath + 'synrefl.sav'
  refl = results
  DEFAULT, NFFT, 32
  N_pair = cece.N_CECE

  df = 1./(NFFT*results.a_over_cs*1e-3/results.interp)
  freq = df*FINDGEN(NFFT/2+1)

  unfilt_xspect = FLTARR(NFFT/2+1)
  unfilt_xspect_err = FLTARR(NFFT/2+1)
  unfilt_coh = FLTARR(NFFT/2+1)
  unfilt_xphase = FLTARR(NFFT/2+1)
  unfilt_xphase_err = FLTARR(NFFT/2+1)
  syn_xspect = FLTARR(NFFT/2+1)
  syn_xspect_err = FLTARR(NFFT/2+1)
  syn_coh = FLTARR(NFFT/2+1)
  syn_xphase = FLTARR(NFFT/2+1)
  syn_xphase_err = FLTARR(NFFT/2+1)

  FOR i_tf=0, results.n_tor_frac-1 DO FOR ir = 0, 2*N_pair-1 DO BEGIN
      unfilt_xspect += CALC_RFSPECT1D_PHASE_ONE(refl.gyro_ne[ir,i_tf,*], $
                                                cece.gyro_te[ir,i_tf,*], $
                                                NFFT=NFFT, ERR=err, $
                                                COHERENCY=coh, PHASE=ph, $
                                                APHAS_ERR=pherr)   
      unfilt_xspect_err += err
      unfilt_coh += coh
      unfilt_xphase += ph
      unfilt_xphase_err += pherr
          
      syn_xspect += CALC_RFSPECT1D_PHASE_ONE(refl.refl[ir,i_tf,*], $
                                             cece.cece[ir,i_tf,*], $
                                             NFFT=NFFT, ERR=err, $
                                             COHERENCY=coh, PHASE=ph, $
                                             APHAS_ERR=pherr)   
      syn_xspect_err += err
      syn_coh += coh
      syn_xphase += ph
      syn_xphase_err += pherr
  ENDFOR
  NM = results.N_tor_frac*2*N_pair
  unfilt_xspect /= NM
  unfilt_xspect_err /= NM
  unfilt_coh /= NM
  unfilt_xphase /= NM
  unfilt_xphase_err /= NM
  unfilt_coh_err = FLTARR(NFFT/2+1)
  syn_xspect /= NM
  syn_xspect_err /= NM
  syn_coh /= NM
  syn_xphase /= NM
  syn_xphase_err /= NM
  syn_coh_err = FLTARR(NFFT/2+1)

  df = 1000./(NFFT*results.a_over_cs/results.interp) ;df in kHz, assume data saved every a/cs
  freq = df*FINDGEN(NFFT/2+1)

  !P.MULTI=[0,0,3]
  PLOT, freq, unfilt_xspect, YTITLE = 'xspectra', CHARSIZE=3, /XS
  ERRPLOT, freq, unfilt_xspect-unfilt_xspect_err, unfilt_xspect + unfilt_xspect_err
  OPLOT, freq, syn_xspect, COLOR=100
  ERRPLOT, freq, syn_xspect-syn_xspect_err, syn_xspect + syn_xspect_err, COLOR=100
  
  PLOT, freq, unfilt_coh, YTITLE = 'coherency', CHARSIZE=3, YRANGE=[0,1], /YS, /XS
;  ERRPLOT, freq, unfilt_coh-unfilt_coh_err, unfilt_coh + unfilt_coh_err
  OPLOT, freq, syn_coh, COLOR=100
;  ERRPLOT, freq, syn_coh-syn_coh_err, syn_coh + syn_coh_err, COLOR=100

  PLOT, freq, unfilt_xphase, YTITLE = 'xphase (rad)', CHARSIZE=3, YRANGE=[-!PI,!PI], /YS,$
        XTITLE = 'frequency (kHZ)', /XS
  OPLOT, freq, FLTARR(NFFT/2+1), LINESTYLE=2
  ERRPLOT, freq, unfilt_xphase-unfilt_xphase_err, unfilt_xphase + unfilt_xphase_err
  OPLOT, freq, syn_xphase, COLOR=100
  ERRPLOT, freq, syn_xphase-syn_xphase_err, syn_xphase + syn_xphase_err, COLOR=100
  !P.MULTI=0

  DEFAULT, fmin, 0.
  DEFAULT, fmax, MAX(freq)
  mind = MIN(ABS(freq-fmin),idx1)
  mind = MIN(ABS(freq-fmax),idx2)

  PRINT, 'df = ', freq[1]
  PRINT, 'fmin = ', freq[idx1]
  PRINT, 'fmax = ', freq[idx2]
  avg_unfilt_xphase = MEAN(unfilt_xphase[idx1:idx2])
  avg_unfilt_xphase_err = MEAN(unfilt_xphase_err[idx1:idx2])
  avg_syn_xphase = MEAN(syn_xphase[idx1:idx2])
  avg_syn_xphase_err = MEAN(syn_xphase_err[idx1:idx2])

  PRINT, 'average unfiltered cross-phase: ', avg_unfilt_xphase ,'+/-', avg_unfilt_xphase_err
  PRINT, 'average synthetic cross-phase: ', avg_syn_xphase ,'+/-', avg_syn_xphase_err
  result={$
   freq:freq,$
   unfilt_xspect:unfilt_xspect,$
   unfilt_xspect_err:unfilt_xspect_err,$
   unfilt_coh:unfilt_coh,$
   unfilt_xphase:unfilt_xphase,$
   unfilt_xphase_err:unfilt_xphase_err,$
   syn_xspect:syn_xspect,$
   syn_xspect_err:syn_xspect_err,$
   syn_coh:syn_coh,$
   syn_xphase:syn_xphase,$
   syn_xphase_err:syn_xphase_err,$
   avg_unfilt_xphase:avg_unfilt_xphase,$
   avg_unfilt_xphase_err:avg_unfilt_xphase_err,$
   avg_syn_xphase:avg_syn_xphase,$
   avg_syn_xphase_err:avg_syn_xphase_err }

  IF N_ELEMENTS(file) NE 0 THEN BEGIN
      OPENW, 1, file
      PRINTF, 1, '#freq unfilt   unfilt_err    syn    syn_err'
      FOR ii=0,NFFT/2 DO PRINTF, 1, freq[ii], $
        unfilt_xphase[ii], unfilt_xphase_err[ii], $
        syn_xphase[ii], syn_xphase_err[ii]
      CLOSE, 1
  ENDIF
END ;analyze_syncece
