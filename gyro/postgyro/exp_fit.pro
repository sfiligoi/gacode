PRO absexp_funct, X, A, F, pder

  F = EXP(-ABS(X)/A)
  IF N_PARAMS() GE 4 THEN pder = (ABS(X)/A^2)*EXP(-ABS(X)/A)
END ;exp_funct

PRO exp_funct, X, A, F, pder

  F = EXP(-X/A)
  IF N_PARAMS() GE 4 THEN pder = (X/A^2)*EXP(-X/A)
END ;exp_funct

FUNCTION exp_fit, x, y, tau_init, ABSX = absx
;
; C. Holland, UCSD
; v1.00: 2/9/2007
;
; fits y(x) with exp(-x/tau), returns tau.  tau_init is an inital guess for
; tau.  No weighting is used.
;
; v2: 7/20/2007: absx flag to fix exp(-|x|/tau)
;
  IF N_ELEMENTS(tau_init) EQ 0 THEN tau = 1. ELSE tau = tau_init

  IF KEYWORD_SET(absx) THEN $
    fit = CURVEFIT(X, Y, weights, tau, FUNCTION = 'absexp_funct') $
  ELSE fit = CURVEFIT(X, Y, weights, tau, FUNCTION = 'exp_funct')
  RETURN, tau
END ;exp_fit
