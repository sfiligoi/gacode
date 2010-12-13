FUNCTION calc_fspect1d, x, y, NFFT=NFFT, FREQUENCY = frequency, $
  ERR = err, NO_HANNING = no_hanning, KEEP_MEAN=keep_mean
;
; C. Holland, 3/28/07
; UCSD
;
; Calculates *MAGNITUDE* of cross-spectrum |S_xy(f)| = |<X(-f)Y(f)>|.
; If nothings is specified for y, defaults to autopower spectrum w/
; y=x.  Assumes x, y are *COMPLEX* 1D arrays of same length, so
; calculates both positive and negative frequencies, returned in
; FREQUENCY keyword (under the assumption of DT = 1)
;
; 10/28/2007: added keep_mean keyword
;
  NT = N_ELEMENTS(x)
  IF N_ELEMENTS(y) EQ 0 THEN y = x ELSE IF (N_ELEMENTS(y) NE NT) THEN BEGIN
      HELP, x, y
      MESSAGE, 'x and y must have same length! Returning 0', /INFO
      RETURN, 0
  ENDIF

  NM = NT/NFFT

  IF KEYWORD_SET(no_hanning) THEN hanwin = FLTARR(NFFT) + 1 $
  ELSE hanwin = HANNING(NFFT)
  hanW = NFFT/TOTAL(hanwin^2)

  CrossSpect = COMPLEXARR(NFFT)
  Sig1Spect = FLTARR(NFFT)
  Sig2Spect = FLTARR(NFFT)
  FOR im = 0L, NM-1L DO BEGIN
	xhat = x[im*NFFT:(im+1)*NFFT-1L]
	yhat = y[im*NFFT:(im+1)*NFFT-1L]
        IF KEYWORD_SET(keep_mean) THEN BEGIN
            xhat = FFT(xhat*hanwin, /INVERSE)/NFFT
            yhat = FFT(yhat*hanwin, /INVERSE)/NFFT
        ENDIF ELSE BEGIN
            xhat = FFT((xhat - MEAN(xhat))*hanwin, /INVERSE)/NFFT
            yhat = FFT((yhat - MEAN(yhat))*hanwin, /INVERSE)/NFFT
        ENDELSE
      	CrossSpect += CONJ(xhat)*yhat
	Sig1Spect += ABS(xhat)^2
	Sig2Spect += ABS(yhat)^2
  ENDFOR
  coherency = SHIFT(ABS(CrossSpect)/(SQRT(Sig1Spect)*SQRT(Sig2Spect)),NFFT/2)

  CrossSpect = SHIFT(ABS(CrossSpect)*hanW/NM,NFFT/2)
  Sig1Spect = SHIFT(Sig1Spect*hanW/NM,NFFT/2)
  Sig2Spect = SHIFT(Sig2Spect*hanW/NM,NFFT/2)
  err = SQRT(Sig1Spect*Sig2Spect*(1. - coherency^2)/NM)
  df = 2*!PI/NFFT
  frequency = df*(FINDGEN(NFFT)-NFFT/2)
  
  RETURN, crossspect
END ;calc_fspect1d

PRO testspect, sf

  DEFAULT, sf, 1
  NT=16
  NM=100
  t = FINDGEN(NT*NM)
  y = sin(2*!PI*4*t/NT)

  yspect = CALC_FSPECT1D(y,NFFT=16*sf,FREQUENCY=frequency)
  yspect2 = CALC_FSPECT1D(y,NFFT=16*sf*2,FREQUENCY=frequency2)
  plot, frequency2, yspect2, psym=-4 ;, yrange=[0,0.25], /xs
  oplot, frequency, yspect, color=100, psym=-4

  print, total(yspect), total(yspect2), total(y^2)/(NM*NT)
END ;testspect

PRO testspect2, NFFT

  DEFAULT, NFFT, 512
  RESTORE, '~/research/BES/BES_DATA/118855_n1-4filt.sav'

  S11 = CALC_FSPECT1D(n1,n1,NFFT=NFFT,FREQ=freq,ERR=err)
  S12 = CALC_FSPECT1D(n1,n2,NFFT=NFFT,ERR=err2)
  freq /= 2*!PI/1E3
  PLOT, freq, S11,XRANGE=[20,250]
  ERRPLOT, freq, S11-err,S11+err
  OPLOT, freq, S12, COLOR=100
  ERRPLOT, freq, S12-err2,S12+err2, COLOR=100

END ;testspec2
