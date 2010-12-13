PRO gauss2_funct_nt1, X, A, F, pder

  F = EXP(-(X/A)^2)
  IF N_PARAMS() GE 4 THEN pder = 2*(X^2/A^3)*EXP(-(X/A)^2)
END ;gauss2_funct_nt1

PRO gauss2_funct_nt2, X, A, F, pder

  F = A[0]*EXP(-(X/A[1])^2)
  dFdA0 = EXP(-(X/A[1])^2)
  dFdA1 = (X^2/A[1]^3)*F
  IF N_PARAMS() GE 4 THEN pder = [[dFdA0],[dFdA1]]

END ;gauss2_funct_nt2

FUNCTION gauss_fit2, x, y, a_init, WEIGHTS=weights, SIGMA=sigma, $
                     CHISQ=chisq, NTERMS=nterms
;
; C. Holland, UCSD
; v1.00: 10/28/2007
;
; IF nterms=1 (default):
; fits y(x) with exp(-(x/a)^2), returns a.  a_init is an inital guess for
; tau.
; IF nterms=2:
; fits y with a[0]*exp(-(x/a[1])^2)
; v2.0: 3/17/2008- added sigma keyword
; v3.0: 5/1/2008- added weights and chisq keywords, and nterms
;

  DEFAULT, nterms, 1
  IF (nterms EQ 1) THEN BEGIN
      IF N_ELEMENTS(a_init) NE 0 THEN a = a_init ELSE a = 1.
      fit = CURVEFIT(X, Y, weights, a, sigma, $
                     FUNCTION_NAME = 'gauss2_funct_nt1', $
                     CHISQ=chisq)
  ENDIF ELSE BEGIN
      IF N_ELEMENTS(a_init) NE 0 THEN a = a_init ELSE a = [1.,1.]
      fit = CURVEFIT(X, Y, weights, a, sigma, $
                     FUNCTION_NAME = 'gauss2_funct_nt2', $
                     CHISQ=chisq)
  ENDELSE
  RETURN, a
END ;gauss_fit2
