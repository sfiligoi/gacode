PRO GYRO_GENERATE_RZCOORDS, profile_data, n_theta, Rmaj, Z
;
; C. Holland, UCSD
;
; Takes rmin and theta (= rho and u in M-L), transforms into 2D arrays
; Rmaj and Z
;
; v1.0: 11/22/06.  Only simple geometry, need to add in rest
; v2.0: 2/16/07.  Added profile_data option for shaping; not passing
;                 this leads to circular flux surfaces.  Defaults for
;                 rmin and theta should be rmin = data.r and 
;                 theta = -!PI + 2*!PI*FINDGEN(n_theta)/ntheta,
;                 specified before calling.  I left it this way for
;                 backwards compatibility, and to keep the routine general
; v3.0: 4/8/10.  Updated for compatibility with current Miller
; forumulation.  General equilibirum shapes not yet supported 
;

  theta = -!PI + 2*!PI*FINDGEN(n_theta)/n_theta
  Ir = FLTARR(profile_data.n_r) + 1
  Ith = FLTARR(n_theta) + 1.

  r2d = profile_data.r # Ith
  theta2d = Ir # theta  
  xd = ASIN(profile_data.delta # Ith)

  Rmaj = (profile_data.R0 # Ith) + r2d*COS(theta2d + xd*SIN(theta2d))
  Z = (profile_data.zmag # Ith) + ((profile_data.kappa # Ith)*r2d*$
       SIN(theta2d + (profile_data.zeta # Ith)*sin(2*theta2d)))
END ;GYRO_GENERATE_RZCOORDS
