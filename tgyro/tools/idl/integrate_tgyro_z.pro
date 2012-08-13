FUNCTION integrate_tgyro_z, prof0, z, rmin, rmin_exp
;
; C. Holland, UCSD
;
; given inv. scale length profile z = -d ln f/d (rmin/a) on
; coarse radial grid rmin, generates f profile on fine rmin_exp grid
; with BC that f(0) = prof0 (does not match at "pivot" point)
;
; sets f=0 at rmin_exp > MAX(rmin)
;
  n_exp = N_ELEMENTS(rmin_exp)
  z_exp = INTERPOL(z, rmin, rmin_exp)
  prof_exp = FLTARR(n_exp)
  ir_max = MAX(WHERE(rmin_exp LE MAX(rmin)))

  FOR ir = 1, ir_max DO $
           prof_exp[ir] = INT_TABULATED(rmin_exp[0:ir],z_exp[0:ir])
  prof_exp = prof0*EXP(-prof_exp)

  prof_exp[ir_max+1:n_exp-1] = 0.

  RETURN, prof_exp  
END ;integrate_tgyro_z

