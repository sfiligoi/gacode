FUNCTION calc_xphase1d, x, y, NFFT=NFFT, FREQUENCY = frequency, $
  ERR = err, NO_HANNING = no_hanning
;
; C. Holland, 3/28/07
; UCSD
;
; Calculates coherency xphase_xy(f) = ATAN(S_xy(f))
; If nothings is specified for y, defaults to y=x.  
; Assumes x, y are *COMPLEX* 1D arrays of same length, so
; calculates both positive and negative frequencies, returned in
; FREQUENCY keyword (under the assumption of DT = 1)
;
; does not assume x and y are real; returns both positive and negative
; frequency values
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
	xhat = FFT((xhat - MEAN(xhat))*hanwin, /INVERSE)/NFFT
	yhat = y[im*NFFT:(im+1)*NFFT-1L]
	yhat = FFT((yhat - MEAN(yhat))*hanwin, /INVERSE)/NFFT
      	CrossSpect += CONJ(xhat)*yhat
	Sig1Spect += ABS(xhat)^2
	Sig2Spect += ABS(yhat)^2
  ENDFOR
  coherency = SHIFT(ABS(CrossSpect)/(SQRT(Sig1Spect)*SQRT(Sig2Spect)),NFFT/2)
  xphase = SHIFT(ATAN(CrossSpect, /PHASE),NFFT/2)

  err = SQRT(1. - coherency^2)/(SQRT(2*NM)*coherency)*!PI
  df = 2*!PI/NFFT
  frequency = df*(FINDGEN(NFFT)-NFFT/2)
  
  RETURN, xphase
END ;calc_xphase1d

PRO testxphase, NFFT

  DEFAULT, NFFT, 512
  RESTORE, '~/research/BES/BES_DATA/118855_n1-4filt.sav'

  xphase = CALC_XPHASE1D(n1,n2,NFFT=NFFT,FREQUENCY=freq,ERR=err)
  freq /= 2*!PI/1E3
  PLOT, freq, xphase,XRANGE=[0,250]
  ERRPLOT, freq, xphase-err,xphase+err

END ;testxphase
