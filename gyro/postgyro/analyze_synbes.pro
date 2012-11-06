PRO analyze_synbes, simdir, NFFT=NFFT, FMIN=fmin, FMAX=fmax, $
                     LOCAL=local
; C. Holland, UCSD
;;
; cacluates statistics for make_synbes_array output, assumign
; output is saved in array named 'results' in file synbes.sav in
; simdir
;
; by default, assumes working on global simulation.  Set local flag to
; average over radius as well

  dirpath = GETENV('GYRO_DIR') + '/sim/' + simdir + '/'
  RESTORE, dirpath + 'synbes.sav'
  DEFAULT, NFFT, 32
  NR_BES = results.NR_BESchannels

  HELP, results, /STRUCT

  df = 1./(NFFT*results.a_over_cs*1e-3/results.interp)
  freq = df*FINDGEN(NFFT/2+1)

  gyro_ne_xspect = FLTARR(NR_BES, NFFT/2+1)
  gyro_ne_xspect_err = FLTARR(NR_BES, NFFT/2+1)
  syn_ne_xspect = FLTARR(NR_BES, NFFT/2+1)
  syn_ne_xspect_err = FLTARR(NR_BES, NFFT/2+1)

  FOR i_tf=0, results.n_tor_frac-1 DO FOR iz=0,results.NZ_BESchannels-2 DO $
    FOR ir = 0, NR_BES-1 DO BEGIN

      gyro_ne_xspect[ir,*] += CALC_RFSPECT1D(results.gyro_ne[ir,iz,i_tf,*], $
                                       results.gyro_ne[ir,iz+1,i_tf,*], $
                                       NFFT=NFFT, ERR=err)
      gyro_ne_xspect_err[ir,*] += err^2

      syn_ne_xspect[ir,*] += CALC_RFSPECT1D(results.syn_ne[ir,iz,i_tf,*], $
                                       results.syn_ne[ir,iz+1,i_tf,*], $
                                       NFFT=NFFT, ERR=err)
      syn_ne_xspect_err[ir,*] += err^2
  ENDFOR
  NM = results.N_tor_frac*(results.NZ_BESchannels-1)
  gyro_ne_xspect /= NM*df
  gyro_ne_xspect_err = SQRT(gyro_ne_xspect_err)/(NM*df)
  syn_ne_xspect /= NM*df
  syn_ne_xspect_err = SQRT(syn_ne_xspect_err)/(NM*df)

  PRINT, 'NM = ', NM

  !P.MULTI = [0,0,2]
  PLOT, freq, gyro_ne_xspect[0,*], /XS, YTITLE = '<|n_e(f)|**2> (1/kHz)',charsize=2
  ERRPLOT, freq, gyro_ne_xspect[0,*]-gyro_ne_xspect_err[0,*], $
           gyro_ne_xspect[0,*]+gyro_ne_xspect_err[0,*]
  OPLOT, freq, gyro_ne_xspect[NR_BES-1,*], LINESTYLE=2
  ERRPLOT, freq, gyro_ne_xspect[NR_BES-1,*]-gyro_ne_xspect_err[NR_BES-1,*], $
           gyro_ne_xspect[NR_BES-1,*]+gyro_ne_xspect_err[NR_BES-1,*]

  OPLOT, freq, syn_ne_xspect[0,*], COLOR=100
  ERRPLOT, freq, syn_ne_xspect[0,*] - syn_ne_xspect_err[0,*], $
           syn_ne_xspect[0,*] + syn_ne_xspect_err[0,*], COLOR=100
  OPLOT, freq, syn_ne_xspect[NR_BES-1,*], LINESTYLE=2, COLOR=100
  ERRPLOT, freq, syn_ne_xspect[NR_BES-1,*]-syn_ne_xspect_err[NR_BES-1,*], $
           syn_ne_xspect[NR_BES-1,*]+syn_ne_xspect_err[NR_BES-1,*], COLOR=100

  DEFAULT, fmin, 0.
  DEFAULT, fmax, MAX(freq)
  mind = MIN(ABS(freq-fmin),idx1)
  mind = MIN(ABS(freq-fmax),idx2)

  PRINT, 'df = ', freq[1]
  PRINT, 'fmin = ', freq[idx1]
  PRINT, 'fmax = ', freq[idx2]

  IF KEYWORD_SET(local) THEN BEGIN
      gyro_ne_xspect = TOTAL(gyro_ne_xspect,1)/NR_BES
      gyro_ne_xspect_err = SQRT(TOTAL(gyro_ne_xspect_err^2,1)/NR_BES)
      syn_ne_xspect = TOTAL(syn_ne_xspect,1)/NR_BES
      syn_ne_xspect_err = SQRT(TOTAL(syn_ne_xspect_err^2,1)/NR_BES)
      
      PLOT, freq, gyro_ne_xspect, xtitle = 'freq (kHz)'
      ERRPLOT, freq, gyro_ne_xspect-gyro_ne_xspect_err, $
               gyro_ne_xspect+gyro_ne_xspect_err
      OPLOT, freq, syn_ne_xspect,color=100
      ERRPLOT, freq, syn_ne_xspect-syn_ne_xspect_err,$
               syn_ne_xspect+syn_ne_xspect_err, color=100

      gyro_ne_rms = SQRT(INT_TABULATED(freq[idx1:idx2],gyro_ne_xspect[idx1:idx2]))
      syn_ne_rms = SQRT(INT_TABULATED(freq[idx1:idx2],syn_ne_xspect[idx1:idx2]))
      gyro_ne_rms_err = 0.5*(INT_TABULATED(freq[idx1:idx2],$
                                             gyro_ne_xspect_err[idx1:idx2]))/$
		    	gyro_ne_rms
      syn_ne_rms_err = 0.5*(INT_TABULATED(freq[idx1:idx2],$
                                              syn_ne_xspect_err[idx1:idx2]))/$
			syn_ne_rms
      PRINT, 'delta S_syn/S_syn: ', INT_TABULATED(freq[idx1:idx2],$
                                                  syn_ne_xspect_err[idx1:idx2])/$
             INT_TABULATED(freq[idx1:idx2],syn_ne_xspect[idx1:idx2])
      PRINT, 'RMS unfiltered n_e (%): ', gyro_ne_rms*100, $
             '   +/- ', gyro_ne_rms_err*100
      PRINT, 'RMS synthetic n_e (%): ', syn_ne_rms*100, $
             '   +/- ', syn_ne_rms_err*100
      
  ENDIF ELSE BEGIN
      
      data = GET_GYRO_DATA(simdir)

      gyro_ne_rms = FLTARR(NR_BES)
      syn_ne_rms = FLTARR(NR_BES)
      gyro_ne_rms_err = FLTARR(NR_BES)
      syn_ne_rms_err = FLTARR(NR_BES)
      
      FOR ir=0,NR_BES-1 DO BEGIN
          gyro_ne_rms[ir]=SQRT(INT_TABULATED(freq[idx1:idx2],gyro_ne_xspect[ir,idx1:idx2]))
          syn_ne_rms[ir]=SQRT(INT_TABULATED(freq[idx1:idx2],syn_ne_xspect[ir,idx1:idx2]))
          gyro_ne_rms_err[ir]=0.5*(INT_TABULATED(freq[idx1:idx2],$
                                                 gyro_ne_xspect_err[ir,idx1:idx2]))/$
                              gyro_ne_rms[ir]
          syn_ne_rms_err[ir] = 0.5*(INT_TABULATED(freq[idx1:idx2],$
                                                  syn_ne_xspect_err[ir,idx1:idx2]))/$
                               syn_ne_rms[ir]
          PRINT, 'r_min/a: ', results.syn_bes_xloc[ir,0]-data.R0[data.n_r/2]
          print, 'delta S_syn / S_syn: ', int_Tabulated(freq[idx1:idx2],$
                                                        syn_ne_xspect_err[ir,idx1:idx2])/$
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


      xloc = results.syn_bes_xloc[*,0] - data.R0[data.n_r/2]
      PLOT, xloc, gyro_ne_rms, psym=4
      ERRPLOT,xloc, gyro_ne_rms-gyro_ne_rms_err,$
              gyro_ne_rms+gyro_ne_rms_err
      
      OPLOT, xloc, syn_ne_rms, PSYM=2, color=100
      ERRPLOT, xloc, syn_ne_rms-syn_ne_rms_err, syn_ne_rms+syn_ne_rms_err, $
               color=100
     
      oplot, xloc, syn_ne_rms*sqrt(1. + 2*syn_ne_rms_err/syn_ne_rms), color=100,linestyle=2
      oplot, xloc, syn_ne_rms*sqrt(1. - 2*syn_ne_rms_err/syn_ne_rms), color=100,linestyle=2
	

	PRINT, 'rmin/a', xloc
	PRINT, 'synBES (%): ', syn_ne_rms
  ENDELSE
  !P.MULTI=0


END ;analyze_syncece
