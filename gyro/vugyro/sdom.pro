FUNCTION SDOM, x
;
; Calculates standard deviation of mean = SQRT(var(x)/N), 
; where N = N_ELEMENTS(x)
;

  RETURN, SQRT(VARIANCE(x)/N_ELEMENTS(x))
END ;SDOM
