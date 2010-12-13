FUNCTION calc_corr1d, x, y, NFFT=NFFT, LAG = lag, NO_HANNING = no_hanning
;
; C. Holland, 12/21/06
; UCSD
;
; Calculates *MAGNITUDE* of cross-correlation 
; C_xy(tau) = |<CONJ(X(t))Y(t+tau)>|/SQRT(<|X(t)|^2><|Y(t)|^2>).
; If nothing is specified for y, defaults to autocorrelation w/
; y=x.  Assumes x, y are *COMPLEX* 1D arrays of same length.
;

  IF N_ELEMENTS(y) EQ 0 THEN y = x
  NT = N_ELEMENTS(x)
  IF N_ELEMENTS(NFFT) EQ 0 THEN NFFT = NT
  IF N_ELEMENTS(NFFT) GT NT THEN NFFT = NT
  lag = FINDGEN(NFFT) - NFFT/2

  NM = NT/NFFT
  hanwin = HANNING(NFFT)
  IF KEYWORD_SET(no_hanning) THEN hanwin = FLTARR(NFFT) + 1
  Sxy = COMPLEXARR(NFFT)
  Px = FLTARR(NFFT)
  Py = FLTARR(NFFT)
  FOR im = 0L, NM-1 DO BEGIN
      idx = im*NFFT + LINDGEN(NFFT)
      psi_x = FFT((x[idx] - MEAN(x[idx]))*hanwin, /INVERSE)
      psi_y = FFT((y[idx] - MEAN(y[idx]))*hanwin, /INVERSE)
      Sxy += CONJ(psi_x)*psi_y
      Px += ABS(psi_x)^2
      Py += ABS(psi_y)^2
  ENDFOR
  Rxy = FLOAT(FFT(Sxy))
  Rx = FLOAT(FFT(Px))
  Ry = FLOAT(FFT(Py))

  Cxy = SHIFT(Rxy/SQRT(Rx[0]*Ry[0]), NFFT/2)
;  Rxy = SHIFT(FLOAT(FFT(Sxy)), NFFT/2)
;  Rx = FLOAT(FFT(Px))
;  Ry = FLOAT(FFT(Py))
;  Cxy = Rxy/SQRT(Rx[0]*Ry[0])
  RETURN, Cxy
END ;calc_corr1d

PRO testcorr, sf

  DEFAULT, sf, 1
  NT=16
  NM=100
  t = FINDGEN(NT*NM)
  y = sin(2*!PI*4*t/NT)

  ycorr = CALC_CORR1D(y,NFFT=16*sf,LAG=lag)
  ycorr2 = CALC_CORR1D(y,NFFT=16*sf*2,LAG=lag2)
  PLOT,  lag2, ycorr2, /XS, PSYM=-4
  OPLOT, lag, ycorr, COLOR=100, PSYM=-2
  OPLOT, lag2, A_CORRELATE(y,lag2), PSYM=-4, COLOR=50
;  env2 = ABS(COMPLEX(ycorr2,HILBERT(ycorr2)))
;  env2 = ABS(ycorr2 + COMPLEX(0.,1.)*HILBERT(ycorr2))
;  PLOT, lag2, env2, COLOR=50

END ;testcorr

PRO testcorr2, NFFT

  DEFAULT, NFFT, 512
  RESTORE, '~/research/BES/BES_DATA/118855_n1-4filt.sav'

  NM = N_ELEMENTS(n1)/NFFT
  dt = 1E-3

  !P.MULTI = [0,0,2,0,0]

  c1 = CALC_CORR1D(n1, n1, NFFT=NFFT,lag=lag)
  c2 = CALC_CORR1D(n1, n1, NFFT=NFFT*2,lag=lag2)
  c3 = CALC_CORR1D(n1, n1, NFFT=NFFT/2,lag=lag3)

  PLOT, lag*dt, c1, $ 
        /XS, XRANGE=[-0.1,0.1], $ 
        /YS, YRANGE=[-1,1], $
        TITLE = 'NFFT = ' + NUMTOSTRING(NFFT) +', NM= ' + NUMTOSTRING(NM), $
        YTITLE = 'C_11(tau)'
  OPLOT, lag2*dt, c2, COLOR=50
  OPLOT, lag3*dt, c3, COLOR=100

  OPLOT, lag*dt, ENVELOPE(c1), LINESTYLE=2
  OPLOT, lag2*dt, ENVELOPE(c2), COLOR=50, LINESTYLE=2
  OPLOT, lag3*dt, ENVELOPE(c3), COLOR=100, LINESTYLE=2

  c1 = CALC_CORR1D(n1, n4, NFFT=NFFT,lag=lag)
  c2 = CALC_CORR1D(n1, n4, NFFT=NFFT*2,lag=lag2)
  c3 = CALC_CORR1D(n1, n4, NFFT=NFFT/2,lag=lag3)

  PLOT, lag*dt, c1, $ 
        /XS, XRANGE=[-0.1,0.1], $ 
        /YS, YRANGE=[-1,1], $
       XTITLE = 'tau (ms)', YTITLE = 'C_14(tau)
        
  OPLOT, lag2*dt, c2, COLOR=50
  OPLOT, lag3*dt, c3, COLOR=100

  OPLOT, lag*dt, ENVELOPE(c1), LINESTYLE=2
  OPLOT, lag2*dt, ENVELOPE(c2), COLOR=50, LINESTYLE=2
  OPLOT, lag3*dt, ENVELOPE(c3), COLOR=100, LINESTYLE=2

  !P.MULTI = 0
END ;testcorr2

