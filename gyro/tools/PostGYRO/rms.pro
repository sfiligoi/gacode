FUNCTION RMS, X
;
; C. Holland, UCSD
; v1.00: 2/6/2008
;
; macro to calculate RMS value of array X

  NX = N_ELEMENTS(X)

  rms = SQRT(TOTAL(X^2)/NX)

  RETURN, rms
END ;RMS
