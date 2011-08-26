PRO GYRO_GENERATE_RZCOORDS, data, n_theta, Rmaj, Z
;
; C. Holland, UCSD
;
; Used data structure from GET_GYRO_DATA to transform
; rmin and theta (= rho and u in M-L) into 2D arrays of
; Rmaj and Z
;
; INPUTS:
;        data: data structure from GET_GYRO_DATA
;        n_theta:  # of theta grid points to use in generating Rmaj
;        and Z
;
; OUTPUTS:
;        Rmaj: [data.n_r x n_theta] array of major radius gridpoints,
;        normed to a
;        Z: [data.n_r x n_theta] array of vertical position
;        gridpoints, normed to a
;
; HISTORY:
; v1.0: 11/22/06.  Only simple geometry, need to add in rest
; v2.0: 2/16/07.  Added profile_data option for shaping; not passing
;                 this leads to circular flux surfaces.  Defaults for
;                 rmin and theta should be rmin = data.r and 
;                 theta = -!PI + 2*!PI*FINDGEN(n_theta)/ntheta,
;                 specified before calling.  I left it this way for
;                 backwards compatibility, and to keep the routine general
; v3.0: 4/8/10.  Updated for compatibility with current Miller
; forumulation.  General equilibirum shapes not yet supported 
; v4.0: 8.25.11  Updated for gacode compatibility
;

  theta = -!PI + 2*!PI*FINDGEN(n_theta)/n_theta
  Ir = FLTARR(data.n_r) + 1
  Ith = FLTARR(n_theta) + 1.

  r2d = data.r # Ith
  theta2d = Ir # theta  
  xd = ASIN(data.delta # Ith)

  Rmaj = (data.R0 # Ith) + r2d*COS(theta2d + xd*SIN(theta2d))
  Z = (data.zmag # Ith) + ((data.kappa # Ith)*r2d*$
       SIN(theta2d + (data.zeta # Ith)*sin(2*theta2d)))

END ;GYRO_GENERATE_RZCOORDS
