;
; C. Holland, UCSD
; v1.0: 2-6-2008
;
;
; Calculates magnitude of cross-power Sxy between two signals x and y, 
; assumed to be real, which is properly normalized such that if x=y,
; MEAN(x) = 0, NFFT=NX, and no hanning window is used, then the
; mean-square value of x, 
;
; x_msv = (sum_i=0,NX-1 x_i^2)/NX 
;
; is equal to sum_i=0,NX/2 Sxy_i
;
; Thus if x and y were voltages sampled every second, the units of Sxy
; would be volt**2/sec
;

FUNCTION calc_rfspect1d, x, y, NFFT=NFFT, FREQ=freq, NO_HANNING=no_hanning, $
  ERR=err, COHERENCY=coherency
;
; 3-13-2008: fixed bug in when abs. value of S taken
;
  NX = N_ELEMENTS(x)
  IF N_ELEMENTS(y) EQ 0 THEN y = x
  IF (N_ELEMENTS(y) NE N_ELEMENTS(x)) THEN BEGIN
      PRINT, 'x and y not same size, returning 0'
      freq = FINDGEN(NX)
      RETURN, FLTARR(NX)
  ENDIF

  DEFAULT, NFFT, NX
  freq = FINDGEN(NFFT/2+1)
  S = COMPLEXARR(NFFT)
  Px = FLTARR(NFFT)
  Py = FLTARR(NFFT)
  NM = NX/NFFT

  IF KEYWORD_SET(no_hanning) THEN hanW = FLTARR(NFFT) + 1 $
  ELSE hanW = HANNING(NFFT)
  hw_norm = NFFT/TOTAL(hanW^2)
  FOR im = 0, NM-1 DO BEGIN
      idx = im*NFFT + LINDGEN(NFFT)
      xhat = (x[idx]-MEAN(x[idx]))*hanW
      yhat = (y[idx]-MEAN(y[idx]))*hanW
      xk = FFT(xhat)
      yk = FFT(yhat)
      S += CONJ(xk)*yk
      Px += ABS(xk)^2
      Py += ABS(yk)^2
  ENDFOR
  S = ABS(S[0:NFFT/2])*hw_norm/NM
  S[1:NFFT/2-1] *= 2

  Px = Px[0:NFFT/2]*hw_norm/NM
  Py = Py[0:NFFT/2]*hw_norm/NM
  Px[1:NFFT/2-1] *= 2
  Py[1:NFFT/2-1] *= 2
  coherency = S/SQRT(Px*Py)
  err = SQRT(Px*Py*(1. - coherency^2)/NM)
  
  RETURN, S
END ;calc_rfspect1d
