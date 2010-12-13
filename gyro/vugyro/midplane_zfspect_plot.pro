PRO midplane_zfspect_plot, plot_log, krmax, NFFT_in, plot_shear, SAVEFILE=savefile, $
  ASCIIFILE=asciifile
;
; C. Holland, UCSD
; v1.0: 5/20/2008
;
; Plots the spectrum of a n=0 fluctuation field at the outboard
; midplane, as a function of k_r and omega (plot_shear = 0).  
;
; Setting plot_shear = 1 plots the shear spectrum 
; |d^2 field/dr^2|^2 = k_r^4 |field|^2
;
; krmax specified max kr value to be plotted.  Routine assumes results
; are from a local sim with periodic radial BC's.
;
; NFFT_in specfies # pts used in frequency FFT, set less than
; (it2-it1+1)/2 to average multiple realizations
;

  common GLOBAL
  common PROFILE_SIM_DATA
  common PLOT_VARIABLES
  common MIDPLANE_DATA

  IF ((NFFT_in LE 1) OR (NFFT_in GT (it2-it1+1))) THEN BEGIN
      MESSAGE, 'NFFT too small/big! Using NFFT=' + NUMTOSTRING(it2-it1+1),/INFO
      NFFT = it2-it1+1
  ENDIF ELSE NFFT = NFFT_in
  NM = (it2-it1+1)/NFFT
  hanwin = HANNING(NFFT)
  dt = (t[n_time-1] - t[0])/(n_time-1)
  freq = (2*!PI/(NFFT*dt))*(FINDGEN(NFFT)-NFFT/2)
  freq_title = '!4x!X (c!Ds!N/a)'
  dkr = (2*!PI/(r[n_r-1] - r[0]))
  all_kr = dkr*(FINDGEN(n_r) - n_r/2) ;not sure where/how vugyro calcs kr_rho

  title = tag[i_pwr]+' n = 0 power spectrum'
  IF (plot_shear EQ 1) THEN BEGIN
      pname = 'midplane_zfshearspect_'+tag[i_pwr]
      title = 'k!S!Dr!R!U4!N * ' + tag[i_pwr]+' n = 0 power spectrum'
  ENDIF ELSE Begin
      pname = 'midplane_zfspect_'+tag[i_pwr]
      title = tag[i_pwr]+' n = 0 power spectrum'
  ENDELSE

  spect = FLTARR(n_r,NFFT)
  hw_norm = NFFT/TOTAL(hanwin^2)
  FOR im=0,NM-1 DO BEGIN
      y = REFORM(pwr[i_pwr,*,0,im*NFFT+LINDGEN(NFFT)])
      y -= MEAN(y)
      y *= (FLTARR(n_r)+1)#hanwin
      ypow = ABS(FFT(y))^2
      spect += SHIFT(ypow,n_r/2,NFFT/2)
  ENDFOR
  spect = REVERSE(spect,2)*hw_norm/NM    ;account for diff in sign convention b/w kr,
                                         ;and omega, correct for power loss
                                         ;from hanning window
  freq = REVERSE(freq)
  IF (plot_shear EQ 1) THEN spect *= (all_kr^4) # (FLTARR(NFFT)+1)  

  ;calculate relative low/high frequency power amounts- denote low freq as 
  ;omega < 0.5 cs/a
  all_pow = 0.;
  lowf_pow = 0.;
  pow_1d = FLTARR(NFFT)
  FOR ii = 0, NFFT-1 DO BEGIN
      pow = INT_TABULATED(all_kr,spect[*,ii])
      pow_1d[ii] = pow
      all_pow += pow
      IF (ABS(freq[ii]) LE 0.5) THEN lowf_pow += pow
  ENDFOR

  IF (plot_shear EQ 1) THEN PRINT, 'Plotting spectrum of n=0 shear k_r^4 |'+tag[i_pwr]+'|^2' $
  ELSE PRINT, 'Plotting spectrum of n=0 |' + tag[i_pwr] + '|^2'
  PRINT, 'faction of spectrum at |omega| <= 0.5 cs/a: ', lowf_pow/all_pow
  PRINT, 'faction of spectrum at |omega| > 0.5 cs/a: ', 1. - lowf_pow/all_pow

  ;rotate to make consistent with earlier plots
  xaxis = freq
  xtitle = freq_title
  xrange = [-!PI,!PI]
  yaxis = all_kr[n_r/2:*]*rho_s
  ytitle = kr_rho_string
  yrange = [0, krmax]
  plot_spect = TRANSPOSE(spect[n_r/2:*,*])
  IF (plot_log EQ 1) THEN BEGIN
      plot_spect = ALOG10(plot_spect > MAX(plot_spect)*1e-8) ;round-off limit
      title = 'log!D10!N ' + title
  ENDIF

  ;do plot
  plot_def_new, pname

  if (c_table_max > 0.0) then begin
      clevels = c_table_min+$
                findgen(nlevels)*(c_table_max-c_table_min)/(nlevels-1.0)
  endif
  loadct, CTab, /SILENT
  
  plot_position = [0.15,0.15,0.85,0.85]
  colorbar_position = FLTARR(4)
  colorbar_position[0] = plot_position[2] + $
    0.02 * (plot_position[2] - plot_position[0])
  colorbar_position[1] = plot_position[1]
  colorbar_position[2] = plot_position[2] + $
    0.05 * (plot_position[2] - plot_position[0])
  colorbar_position[3] = plot_position[3]
      
  contour,plot_spect,xaxis,yaxis,$
          nlevels=nlevels,$
          levels=clevels,$
          xtitle=xtitle, $
          xstyle=1,$
          xrange=xrange,$
          yrange=yrange, $
          ytitle=ytitle, $
          ystyle=1,$
          title=title,$
          /fill, $
          position=plot_position

  smin = MIN(plot_spect,MAX=smax)
;  COLORBAR, NCOLORS = !D.TABLE_SIZE, MINRANGE = smin, $
;            MAXRANGE = smax,/VERTICAL, /RIGHT, $
;            POSITION = colorbar_position, FORMAT = '(G10.2)'

  XYOUTS, 0.96, 0.025, "NFFT="+NUMTOSTRING(NFFT)+', NM='+NUMTOSTRING(NM), /NORMAL, $
          ALIGN=1

  plot_finish
  set_line_colors

  ;write .sav file
  IF N_ELEMENTS(savefile) NE 0 THEN BEGIN
      SAVE, all_kr, freq, spect, FILENAME = savefile
      MESSAGE, 'Wrote to ' + savefile, /INFO
  ENDIF

  ;write text file
  IF N_ELEMENTS(asciifile) NE 0 THEN BEGIN
      OPENW, lun, asciifile, /GET_LUN
      PRINTF, lun, n_r, NFFT
      PRINTF, lun, all_kr, freq, spect
      FREE_LUN, lun
      MESSAGE, 'Wrote to ' + asciifile, /INFO
  ENDIF

END ;midplane_zfspect_plot
