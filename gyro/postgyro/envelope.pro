FUNCTION envelope, f
;
; C. Holland, UCSD
; v1.0: 2/9/07
;
; Use Hilbert transform to calculate envelope of function f
;

  cf = COMPLEX(f, HILBERT(f))
  envelope = ABS(cf)

  RETURN, envelope
END ;envelope
