; AE White 07-13-09
; hunt down any differences
;in the SIGN of the phase angle that may be showing up
; Analysis must be checked because measurements of 
; lemo polarity at digitizer showed no switching in hardware.
;
;

;AE WHite 04-03-08
; Added the cross phase output for CECE REFL studies

; C. Holland, UCSD
; v1.0: 2-6-2008
;
;
; Calculates magnitude of cross-power Sxy between two signals x and y, 
; assumed to be real, which is properly normalized such that if x=y,
; MEAN(x) = 0, NFFT=NX, and no hanning windo is used, then the
; mean-square value of x, 
;
; x_msv = (sum_i=0,NX-1 x_i^2)/NX 
;
; is equal to sum_i=0,NX/2 Sxy_i
;
; Thus if x and y were voltages sampled every second, the units of Sxy
; would be volt**2/sec
;


FUNCTION calc_rfspect1d_phase_one, x, y, NFFT=NFFT, FREQ=freq, NO_HANNING=no_hanning, $
  ERR=err, COHERENCY=coherency, PHASE=PHASE, APHAS_ERR=phase_err, ANGLE1=angle1

;FUNCTION calc_rfspect1d_phase_anne, x, y, NFFT=NFFT, FREQ=freq, NO_HANNING=no_hanning, $
 ; ERR=err, COHERENCY=coherency, PHASE=PHASE, ANGLE1=angle1
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
      xk = FFT(xhat, -1)
      yk = FFT(yhat, -1)
      S += CONJ(xk)*yk
     
      Px += ABS(xk)^2
      Py += ABS(yk)^2
   ENDFOR
   ;Bendat and Piersol Random Data 2000 pp 445 Eqn. 11.149
  ;Gxy_k= 2/(N delta t) x [X_k* x Y_k] k= 0...N/2 
  ;Gxy=  S[0:NFFT/2]
  Gxy=S
  angle1 = ATAN(Gxy, /phase)
  phase = ATAN(Gxy[0:NFFT/2], /phase)
  phase = angle1

;angle1 = angle1[0:NFFT/2]*hw_norm/NM
  ;phase = phase[0:NFFT/2]*hw_norm/NM    



  S = ABS(S[0:NFFT/2])*hw_norm/NM
  S[1:NFFT/2-1] *= 2
  
  Px = Px[0:NFFT/2]*hw_norm/NM
  Py = Py[0:NFFT/2]*hw_norm/NM
  Px[1:NFFT/2-1] *= 2
  Py[1:NFFT/2-1] *= 2
  coherency = S/SQRT(Px*Py)
  err = SQRT(Px*Py*(1. - coherency^2)/NM)
   
  phase_err = sqrt(1. - coherency^2)/(coherency * 2 * sqrt(NM))
  
 
  RETURN, S
END ;calc_rfspect1d
