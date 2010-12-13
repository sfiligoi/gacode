FUNCTION gyro_rtheta_transform, field, profile_data, sf, THETA = theta, $
  NU_SILENT = nu_silent, TOR_FRAC=tor_frac
;
; C. Holland, UCSD
;
; v1.0: 2/21/2007: Implements existing VuGYRO algorithm
; v2.0: 2/27/2007: Updated to use nu(theta,r) instead of just
;                  nu = -q*theta.  Assuming it will be included in
;                  profile data in future version.
; v2.1: 4/3/2007: file added in read-in from geometry_arrays.out.  Three
;                 new fields in profile_data: nu_geo, theta_mult, and 
;                 exists_nu_geo.  Use nu geo if exists_nu_geo=1 and 
;                 sf = theta_mult, otherwise use -qtheta.  Warning
;                 sounds in this case, disabled by NU_SILET = 1
;
; v3: 12/21/2007: added in tor_frac keyword.  Allows user to create
;                 at toroidal angle varphi = tor_frac*2*pi/n_dn.
;                 Setting eq =1 should be same as 0.  Useful for
;                 getting extra statistics out of a run (get data at
;                 more locations than just varphi = 0).
;
; takes field = COMPLEXARR(n_theta_plot, n_r,n_n) and maps into
; field_rtheta = FLTARR(n_y=sf*n_theta_plot,n_r)
;
; assume zeta_phys = 0, so have, e.g.
; phi(theta,r) = Re sum_on_n phi_hat(theta,r,n)exp(-2*pi*i*n*q*theta)
; note phi_hat is not periodic in theta
;
; profile_data is result of read_gyro_profile_data; included in result
; from get_gyro_data (as <result>.profile_data)
; SF is interpolation factor; theta is an axis variable corresponding
; to the "generalized" theta angle (i.e. it isn't the physical theta
; in a shaped plasma)
;

  n_r = profile_data.n_r
  n_theta_plot = profile_data.n_theta_plot
  n_n = profile_data.n_n
  n = profile_data.n_0 + profile_data.n_dn*INDGEN(profile_data.n_n)
  C_I = COMPLEX(0,1)

  n_y = n_theta_plot*sf*1L
  field_rtheta = FLTARR(n_y, n_r)
  theta = -!PI + 2*!PI*FINDGEN(n_y)/n_y

  DEFAULT, nu_silent, 0
  IF (profile_data.exists_nu_geo AND (profile_data.theta_mult EQ sf)) THEN $
    nu = profile_data.nu_geo $
  ELSE BEGIN
      nu = -theta # profile_data.q
      IF (nu_silent EQ 0) THEN MESSAGE, "Using nu = -q*theta", /INFO
  ENDELSE
  alpha = nu
  
  ;add in toroidal angle to alpha
  DEFAULT, tor_frac, 0
  alpha += tor_frac*2*!PI/profile_data.n_dn

  Utemp = COMPLEXARR(n_theta_plot+1,n_r)
  Uinterp = COMPLEXARR(n_y+1)
  U = COMPLEXARR(n_y,n_r,n_n)
  phase = COMPLEXARR(n_r)

  ;map from n_theta_plot to n_y
  FOR i_n = 0, n_n-1 DO BEGIN
      Utemp[0:n_theta_plot-1,*] = field[0:n_theta_plot-1,*,i_n]

      ;apply theta BC: z_n(r,,2*pi) = z_n(r,0)exp(I*n*(nu(r,2*pi)-nu(r,0)))
      ;with nu(r,2*pi) - nu(r,0) = -2*pi*q by definition
      phase[*] = EXP(-2*!PI*C_I*n[i_n]*profile_data.q[*])
      Utemp[n_theta_plot,*] = field[0,*,i_n]*phase

      ;interpolate from n_theta_plot to n_y
      FOR i_r = 0, n_r-1 DO BEGIN
          ;Uinterp = INTERPOL(Utemp[*,i_r],n_y+1,/SPLINE)
          Uinterp = INTERPOL(Utemp[*,i_r],n_y+1)
          U[0:n_y-1,i_r,i_n] = Uinterp[0:n_y-1]
      ENDFOR
  ENDFOR

  ;now transform to physical space
  IF (n[0] EQ 0) THEN BEGIN
      field_rtheta[*,*] = FLOAT(U[*,*,0])
      FOR i_n = 1, n_n-1 DO $
                field_rtheta[*,*] += 2*FLOAT(U[*,*,i_n]*EXP(-C_I*n[i_n]*alpha))
  ENDIF ELSE FOR i_n = 0, n_n-1 DO $
                field_rtheta[*,*] += 2*FLOAT(U[*,*,i_n]*EXP(-C_I*n[i_n]*alpha))

  RETURN, field_rtheta

END ;gyro_rtheta_transform
