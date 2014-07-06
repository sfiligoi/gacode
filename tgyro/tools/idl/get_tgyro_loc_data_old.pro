FUNCTION get_tgyro_loc_data_old, simdir, DIRLOC=dirloc, N_ion=N_ion
;
; C. Holland, UCSD
;
; v1.0: Oct. 29, 2008
;
; Reads in input.profiles and various local TGYRO output files,
; returns everything in one large IDL structure.
;
; Usage: To load a simulation, e.g. 'wtest_1', enter:
; IDL> data = get_tgyro_loc_data('wtest_1')
;
; this will return an IDL structure.  To examine the contents of data, enter
; IDL> help, data, /STRUCT
;
; There are two sets of data: arrays of length EXP_N_RHO corresponding to the
; experimental profiles, and arrays of size [NX,N_it+1] corresponding to TGYRO
; predictions vs. radius and iteration number.
;
; details of units etc. can be found at the bottom of this routine where the
; structure is declared.  As much as possible, I've attempted to keep
; the same naming conventions as the TGYRO output files
;
; KEYWORDS
; simdir: string containing name of valid directory in local directory
; dirloc: string with pathname to directory containing simdir
; N_ion: number of dynamic ion species to load.  Defaults to 1.
;
; v1.1: Apr. 15, 2009
; Added support for gamma_p in gradients.out, no_gammap_flag for
; backwards compatibility
;
; v1.2: May 10, 2009
; Added ne_on flag to support evolving n_e profiles.
;
; v1.3: July 20,2009
; Added no_mom flag for backwards compatibility with no-momentum cases
;
; v1.4: Aug 17, 2009
; Revised for compabitiblity with 1.1 version file structure.  No longer
; backwards compatible
;
; v1.4.1: Aug 24, 2009
; Added 3rd ion species from INPUT_profiles.  Removed no_mom, ne_on flags.
;
; 10.18.09: added geometry_extra.out, including direct read of d(vol)/dr
;
; v2.0: July 19, 2010
; added dirloc keyword to allow simulations located in other than $TGYRO_DIR/sim
;
; v3.0: March 3, 2011
; updated to be compatible with new gacode file structure
;
; v3.1: May 2, 2012
; updated to read in expt. data points from files name, ne_data.txt,
; te_data.txt, ti_data.txt in TGYRO dir and pass along- no problem if
;                                                       not available

  IF N_ELEMENTS(simdir) EQ 0 THEN BEGIN
      MESSAGE, 'Need to specify a directory!', /INFO
      RETURN, 0
  ENDIF


;  IF N_ELEMENTS(dirloc) EQ 0 THEN dirpath = simdir + '/' $
;  ELSE dirpath = dirloc + '/' + simdir + '/'
  DEFAULT, dirloc, '.'
  dirpath = dirloc + '/' + simdir + '/'

  ;first read input.profiles to get inital info on fine grid
  s = STRING('#')
  OPENR, 1, dirpath + 'input.profiles', ERROR=err
  IF (err NE 0) THEN BEGIN
      PRINT, "Couldn't open " + dirpath + "input.profiles!  Returning 0"
      RETURN, 0
  ENDIF

  WHILE((STRPOS(s,'#') EQ 0) OR (STRPOS(s,'V') EQ 0) OR (STRPOS(s,'O') EQ 0)) DO BEGIN
      READF, 1, s
  ENDWHILE
  exp_n_rho = FIX(STRMID(s,6,3))
  READF, 1, s
  READF, 1, s
  arr = FLTARR(5,exp_n_rho)

  WHILE (STRPOS(s, 'rho') EQ -1) DO READF, 1, s
  READF, 1, arr
  exp_rho = REFORM(arr[0,*])
  exp_rmin = REFORM(arr[1,*])
  a = exp_rmin[exp_n_rho-1]
  exp_rmin /= a
  exp_Rmaj = REFORM(arr[2,*])/a
  exp_q = REFORM(arr[3,*])
  exp_kappa = REFORM(arr[4,*])

  s = '#'
  WHILE (STRPOS(s, 'delta') EQ -1) DO READF, 1, s
  READF, 1, arr
  exp_delta = REFORM(arr[0,*])
  exp_Te = REFORM(arr[1,*])
  exp_ne = REFORM(arr[2,*])
  exp_Zeff = REFORM(arr[3,*])
  exp_omega0 = REFORM(arr[4,*])

  s = '#'
  WHILE (STRPOS(s, 'pow_e') EQ -1) DO READF, 1, s
  READF, 1, arr
  exp_flow_mom = REFORM(arr[0,*])  ;Nm
  exp_Pe = REFORM(arr[1,*])        ;Mw
  exp_Pi = REFORM(arr[2,*])        ;Mw
  exp_Pexch = REFORM(arr[3,*])     ;Mw
  exp_zeta = REFORM(arr[4,*])

  READF, 1, s
  READF, 1, s
  READF, 1, arr
  exp_flow_beam = REFORM(arr[0,*])  ;units?? MW/keV?
  exp_flow_wall = REFORM(arr[1,*])
  exp_flow = exp_flow_beam + exp_flow_wall
  exp_zmag = REFORM(arr[2,*])   ;m
  exp_ptot = REFORM(arr[3,*])   ;Pa
  ;column 4 unused

  READF, 1, s
  READF, 1, s
  READF, 1, arr
  exp_ni = REFORM(arr[0,*])    ;10**19/m**3
  exp_ni2 = REFORM(arr[1,*])
  exp_ni3 = REFORM(arr[2,*])
  exp_ni4 = REFORM(arr[3,*])
  exp_ni5 = REFORM(arr[4,*])

  READF, 1, s
  READF, 1, s
  READF, 1, arr
  exp_Ti = REFORM(arr[0,*])
  exp_Ti2 = REFORM(arr[1,*])
  exp_Ti3 = REFORM(arr[2,*])
  exp_Ti4 = REFORM(arr[3,*])
  exp_Ti5 = REFORM(arr[4,*])

  READF, 1, s
  READF, 1, s
  READF, 1, arr
  exp_Vtor = REFORM(arr[0,*])
  exp_Vtor2 = REFORM(arr[1,*])
  exp_Vtor3 = REFORM(arr[2,*])
  exp_Vtor4 = REFORM(arr[3,*])
  exp_Vtor5 = REFORM(arr[4,*])

  READF, 1, s
  READF, 1, s
  READF, 1, arr
  exp_Vpol = REFORM(arr[0,*])
  exp_Vpol2 = REFORM(arr[1,*])
  exp_Vpol3 = REFORM(arr[2,*])
  exp_Vpol4 = REFORM(arr[3,*])
  exp_Vpol5 = REFORM(arr[4,*])

  CLOSE, 1

  ;get NX and N_it = # of iterations
;  OPENR, 1, dirpath + 'control.out'
  OPENR, 1, dirpath + 'out.tgyro.control'
  NX = 0  ;force to be integer
  READF, 1, NX
  N_fields = 0 ;# of fields evolved
  READF, 1, N_fields
  N_it = 0
  READF, 1, N_it
  CLOSE, 1
  N_it++  ;add in inital timepoint

  ;read geometry.out- assume doesn't change during iteations
  rmin = FLTARR(NX)      ;rmin/a
  rho = FLTARR(NX)    ;normalized toroidal flux
  q = FLTARR(NX)
  s_hat = FLTARR(NX)  ;(r_min/q)*dq/dr_min
  kappa = FLTARR(NX)
  s_kappa = FLTARR(NX) ;(r_min/kappa)*dkappa/dr_min
  delta = FLTARR(NX)
  s_delta = FLTARR(NX) ;(r_min/delta)*ddelta/dr_min
  shift = FLTARR(NX)
  Rmaj = FLTARR(NX)    ;Rmaj/a
  b_unit = FLTARR(NX)  ;b_unit (T)
  s = ' '
  arr = fltarr(11,NX)
  OPENR, 1, dirpath + 'out.tgyro.geometry.1', ERROR=err
  IF (err NE 0) THEN BEGIN
      PRINT, 'Using old geometry.out'
      OPENR, 1, dirpath + 'geometry.out'
  ENDIF
  READF, 1, s
  READF, 1, s
  READF, 1, arr
  rmin = REFORM(arr[0,*])
  rho = REFORM(arr[1,*])
  q = REFORM(arr[2,*])
  s_hat = REFORM(arr[3,*])
  kappa = REFORM(arr[4,*])
  s_kappa = REFORM(arr[5,*])
  delta = REFORM(arr[6,*])
  s_delta = REFORM(arr[7,*])
  shift = REFORM(arr[8,*])
  Rmaj = REFORM(arr[9,*])
  b_unit = REFORM(arr[10,*])
  CLOSE, 1

  ;read geometry_extra.out
  s = ' '
  OPENR, 1, dirpath + 'out.tgyro.geometry.2', ERROR=err
  arr = FLTARR(9,NX)
  IF (err NE 0) THEN BEGIN
      PRINT, 'Using old geometry_extra.out'
      OPENR, 1, dirpath + 'geometry_extra.out'
      arr = FLTARR(7,NX)
  ENDIF
  READF, 1, s
  READF, 1, s
  READF, 1, arr
  CLOSE, 1
  zmag = REFORM(arr[1,*])
  dzmag = REFORM(arr[2,*])
  zeta = REFORM(arr[3,*])
  s_zeta = REFORM(arr[4,*])
  vol = REFORM(arr[5,*])
  dvoldr = REFORM(arr[6,*])

  ;read chi_e.out
  De_neo = FLTARR(NX,N_it)
  De_tur = FLTARR(NX,N_it)
  chie_neo = FLTARR(NX,N_it)
  chie_tur = FLTARR(NX,N_it)
  s = ' '
  arr = fltarr(5,NX)
;  OPENR, 1, dirpath + 'chi_e.out'
;  FOR ii = 0, N_it-1 DO BEGIN
;      READF, 1, s
;      READF, 1, s
;      READF, 1, arr
;      De_neo[*,ii] = arr[1,*]
;      De_tur[*,ii] = arr[2,*]
;      chie_neo[*,ii] = arr[3,*]
;      chie_tur[*,ii] = arr[4,*]
;  ENDFOR
  CLOSE, 1

  ;read chi_i.out
  Di_neo = FLTARR(NX,N_it)
  Di_tur = FLTARR(NX,N_it)
  chii_neo = FLTARR(NX,N_it)
  chii_tur = FLTARR(NX,N_it)
  s = ' '
  arr = fltarr(5,NX)
;  OPENR, 1, dirpath + 'chi_i.out'
;  FOR ii = 0, N_it-1 DO BEGIN
;      READF, 1, s
;      READF, 1, s
;      READF, 1, arr
;      Di_neo[*,ii] = arr[1,*]
;      Di_tur[*,ii] = arr[2,*]
;      chii_neo[*,ii] = arr[3,*]
;      chii_tur[*,ii] = arr[4,*]
;  ENDFOR
;  CLOSE, 1

  ;read flux_e.out
  Ge_neo = FLTARR(NX,N_it)
  Ge_tur = FLTARR(NX,N_it)
  Qe_neo = FLTARR(NX,N_it)
  Qe_tur = FLTARR(NX,N_it)
  s = ' '
  arr = fltarr(5,NX)
;  OPENR, 1, dirpath + 'flux_e.out'
  OPENR, 1, dirpath + 'out.tgyro.flux_e'
  FOR ii = 0, N_it-1 DO BEGIN
      READF, 1, s
      READF, 1, s
      READF, 1, arr
      Ge_neo[*,ii] = arr[1,*]
      Ge_tur[*,ii] = arr[2,*]
      Qe_neo[*,ii] = arr[3,*]
      Qe_tur[*,ii] = arr[4,*]
  ENDFOR
  CLOSE, 1

  ;read flux_i.out
  Gi_neo = FLTARR(NX,N_it)
  Gi_tur = FLTARR(NX,N_it)
  Qi_neo = FLTARR(NX,N_it)
  Qi_tur = FLTARR(NX,N_it)
  s = ' '
  arr = FLTARR(5,NX)
;  OPENR, 1, dirpath + 'flux_i.out'
  OPENR, 1, dirpath + 'out.tgyro.flux_i'
  FOR ii = 0, N_it-1 DO BEGIN
      READF, 1, s
      READF, 1, s
      READF, 1, arr
      Gi_neo[*,ii] = arr[1,*]
      Gi_tur[*,ii] = arr[2,*]
      Qi_neo[*,ii] = arr[3,*]
      Qi_tur[*,ii] = arr[4,*]
  ENDFOR
  CLOSE, 1

  ;read flux_target.out
  Qi_tot = FLTARR(NX,N_it)
  Qi_target = FLTARR(NX,N_it)
  Qe_tot = FLTARR(NX,N_it)
  Qe_target = FLTARR(NX,N_it)
  Ge_tot = FLTARR(NX,N_it)
  Ge_target = FLTARR(NX, N_it)
  s = ' '
  arr = fltarr(7,NX)
;  OPENR, 1, dirpath + 'flux_target.out'
  OPENR, 1, dirpath + 'out.tgyro.flux_target'
  FOR ii = 0, N_it-1 DO BEGIN
      READF, 1, s
      READF, 1, s
      READF, 1, arr
      Qi_tot[*,ii] = arr[1,*]
      Qi_target[*,ii] = arr[2,*]
      Qe_tot[*,ii] = arr[3,*]
      Qe_target[*,ii] = arr[4,*]
      Ge_tot[*,ii] = arr[5,*]
      Ge_target[*,ii] = arr[6,*]
  ENDFOR
  CLOSE, 1

  ;read mflux_e.out
  Mflux_e_neo = FLTARR(NX,N_it)
  Mflux_e_tur = FLTARR(NX,N_it)
  Expwd_e_tur = FLTARR(NX,N_it) ;turbulent exchange
  s = ' '
  arr = fltarr(4,NX)
;  OPENR, 1, dirpath + 'mflux_e.out'
  OPENR, 1, dirpath + 'out.tgyro.mflux_e'
  FOR ii = 0, N_it-1 DO BEGIN
      READF, 1, s
      READF, 1, s
      READF, 1, arr
      Mflux_e_neo[*,ii] = arr[1,*]
      Mflux_e_tur[*,ii] = arr[2,*]
      Expwd_e_tur[*,ii] = arr[3,*]
  ENDFOR
  CLOSE, 1

  ;read mflux_i.out
  Mflux_i_neo = FLTARR(NX,N_it)
  Mflux_i_tur = FLTARR(NX,N_it)
  Expwd_i_tur = FLTARR(NX,N_it) ;turbulent exchange
  s = ' '
  arr = fltarr(4,NX)
;  OPENR, 1, dirpath + 'mflux_i.out'
  OPENR, 1, dirpath + 'out.tgyro.mflux_i'
  FOR ii = 0, N_it-1 DO BEGIN
      READF, 1, s
      READF, 1, s
      READF, 1, arr
      Mflux_i_neo[*,ii] = arr[1,*]
      Mflux_i_tur[*,ii] = arr[2,*]
      Expwd_i_tur[*,ii] = arr[3,*]
  ENDFOR
  CLOSE, 1

  ;read mflux_target.out
  Mflux_tot = FLTARR(NX,N_it)
  Mflux_target = FLTARR(NX,N_it)
  s = ' '
  arr = fltarr(3,NX)
;  OPENR, 1, dirpath + 'mflux_target.out'
  OPENR, 1, dirpath + 'out.tgyro.mflux_target'
  FOR ii = 0, N_it-1 DO BEGIN
      READF, 1, s
      READF, 1, s
      READF, 1, arr
      Mflux_tot[*,ii] = arr[1,*]
      Mflux_target[*,ii] = arr[2,*]
  ENDFOR
  CLOSE, 1

  ;read gyrobohm.out
  chi_gB = FLTARR(NX,N_it) ;chi_gB in m**/s
  Q_gB = FLTARR(NX,N_it)   ;Q_gB in MW/m**2
  G_gB = FLTARR(NX,N_it)   ;G_gB in 10**19/m**2/s
  Pi_gB = FLTARR(NX,N_it)  ;Pi_gB in J/m**2
  c_s = FLTARR(NX,N_it)    ;c_s in m/s
  s = ' '
  arr=FLTARR(6,NX)

;  OPENR, 1, dirpath + 'gyrobohm.out'
  OPENR, 1, dirpath + 'out.tgyro.gyrobohm'
  FOR ii = 0, N_it-1 DO BEGIN
      READF, 1, s
      READF, 1, s
      READF, 1, arr
      chi_gB[*,ii] = arr[1,*]
      Q_gB[*,ii] = arr[2,*]
      G_gB[*,ii] = arr[3,*]
      Pi_gB[*,ii] = arr[4,*]
      c_s[*,ii] = arr[5,*]
  ENDFOR
  CLOSE, 1

  ;read nu_rho.out
  nu_ii = FLTARR(NX,N_it)       ;(a/cs)/t_ii
  nu_ee = FLTARR(NX,N_it)       ;(a/cs)/t_ee
  nu_estar_inv = FLTARR(NX,N_it)       ;1/((a/cs)nu_e_Star)
  nu_exch = FLTARR(NX,N_it)     ;(a/cs)nu_exch
  rho_i_star = FLTARR(NX,N_it)  ;rho_i/a
  rho_star = FLTARR(NX,N_it)    ;rho_s/a
  s = ' '
  arr = fltarr(8,NX)  ;leave 8th column frac_ae out
  OPENR, 1, dirpath + 'out.tgyro.nu_rho'
  FOR ii = 0, N_it-1 DO BEGIN
      READF, 1, s
      READF, 1, arr
      nu_ii[*,ii] = arr[1,*]
      nu_ee[*,ii] = arr[2,*]
      nu_estar_inv[*,ii] = arr[3,*]
      nu_exch[*,ii] = arr[4,*]
      rho_i_star[*,ii] = arr[5,*]
      rho_star[*,ii] = arr[6,*]
  ENDFOR
  CLOSE, 1

  ;read gradient.out
  a_over_Lni = FLTARR(NX,N_it)
  a_over_Lne = FLTARR(NX,N_it)
  a_over_LTi = FLTARR(NX,N_it)
  a_over_LTe = FLTARR(NX,N_it)
  a_over_Lp = FLTARR(NX,N_it)
  gamma_e = FLTARR(NX,N_it)   ;a*gamma_e/c_s
  gamma_p = FLTARR(NX,N_it)   ;a*gamma_e/c_s
  s = ' '
  arr = fltarr(8,NX)
;  OPENR, 1, dirpath + 'gradient.out'
  OPENR, 1, dirpath + 'out.tgyro.gradient'
  FOR ii = 0, N_it-1 DO BEGIN
      READF, 1, s
      READF, 1, s
      READF, 1, arr
      a_over_Lni[*,ii] = arr[1,*]
      a_over_Lne[*,ii] = arr[2,*]
      a_over_LTi[*,ii] = arr[3,*]
      a_over_LTe[*,ii] = arr[4,*]
      a_over_Lp[*,ii] = arr[5,*]
      gamma_e[*,ii] = arr[6,*]
      gamma_p[*,ii] = arr[7,*]
  ENDFOR
  CLOSE, 1

  ;read profile.out
  n_i = FLTARR(NX,N_it)
  n_e = FLTARR(NX,N_it)
  Ti = FLTARR(NX,N_it)
  Te = FLTARR(NX,N_it)
  betae_unit = FLTARR(NX,N_it)
  M = FLTARR(NX,N_it)
  s = ' '
  arr = fltarr(8,NX)
;  OPENR, 1, dirpath + 'profile.out'
  OPENR, 1, dirpath + 'out.tgyro.profile'
  FOR ii = 0, N_it-1 DO BEGIN
      READF, 1, s
      READF, 1, s
      READF, 1, arr
      n_i[*,ii] = arr[1,*]
      n_e[*,ii] = arr[2,*]
      Ti[*,ii] = arr[3,*]
      Te[*,ii] = arr[4,*]
      betae_unit[*,ii] = arr[6,*]
      M[*,ii] = arr[7,*]
  ENDFOR
  CLOSE, 1

  ;read power.out
  p_alpha = FLTARR(NX,N_it)
  p_brem = FLTARR(NX,N_it)
  p_exch = FLTARR(NX,N_it)
  p_expwd = FLTARR(NX,N_it)
  p_i_aux = FLTARR(NX,N_it)
  p_e_aux = FLTARR(NX,N_it)
  p_i = FLTARR(NX,N_it)
  p_e = FLTARR(NX,N_it)
  s = ' '
  OPENR, 1, dirpath + 'out.tgyro.power', ERROR=err
  IF (err EQ 0 ) THEN BEGIN
;      arr = fltarr(9,NX)
;      FOR ii = 0, N_it-1 DO BEGIN
;          READF, 1, s
;          READF, 1, s
;          READF, 1, arr
;          p_alpha[*,ii] = arr[1,*]
;          p_brem[*,ii] = arr[2,*]
;          p_exch[*,ii] = arr[3,*]
;          p_expwd[*,ii] =arr[4,*]
;          p_i_aux[*,ii] = arr[5,*]
;          p_e_aux[*,ii] = arr[6,*]
;          p_i[*,ii] = arr[7,*]
;          p_e[*,ii] = arr[8,*]
;      ENDFOR
  ENDIF ELSE BEGIN
;      arr = fltarr(8,NX)      
;      PRINT, 'Using old power.out file'
;      OPENR, 1, dirpath + 'power.out'
;      FOR ii = 0, N_it-1 DO BEGIN
;          READF, 1, s
;          READF, 1, s
;          READF, 1, arr
;          p_alpha[*,ii] = arr[1,*]
;          p_brem[*,ii] = arr[2,*]
;          p_exch[*,ii] = arr[3,*]
;          p_i_aux[*,ii] = arr[4,*]
;          p_e_aux[*,ii] = arr[5,*]
;          p_i[*,ii] = arr[6,*]
;          p_e[*,ii] = arr[7,*]
;      ENDFOR
  ENDELSE
  CLOSE, 1

  DEFAULT, N_ion, 1
  Di2_neo = FLTARR(NX,N_it)
  Di2_tur = FLTARR(NX,N_it)
  chii2_neo = FLTARR(NX,N_it)
  chii2_tur = FLTARR(NX,N_it)
  Gi2_neo = FLTARR(NX,N_it)
  Gi2_tur = FLTARR(NX,N_it)
  Qi2_neo = FLTARR(NX,N_it)
  Qi2_tur = FLTARR(NX,N_it)

  IF (N_ion GE 2) THEN BEGIN

     ;read chi_i2.out
     Di2_neo = FLTARR(NX,N_it)
     Di2_tur = FLTARR(NX,N_it)
     chii2_neo = FLTARR(NX,N_it)
     chii2_tur = FLTARR(NX,N_it)
     s = ' '
     arr = fltarr(5,NX)
;     OPENR, 1, dirpath + 'chi_i2.out'
;     FOR ii = 0, N_it-1 DO BEGIN
;        READF, 1, s
;        READF, 1, s
;        READF, 1, arr
;        Di2_neo[*,ii] = arr[1,*]
;        Di2_tur[*,ii] = arr[2,*]
;        chii2_neo[*,ii] = arr[3,*]
;        chii2_tur[*,ii] = arr[4,*]
;     ENDFOR
;     CLOSE, 1

     ;read flux_i2.out
     Gi2_neo = FLTARR(NX,N_it)
     Gi2_tur = FLTARR(NX,N_it)
     Qi2_neo = FLTARR(NX,N_it)
     Qi2_tur = FLTARR(NX,N_it)
     s = ' '
     arr = FLTARR(5,NX)
;     OPENR, 1, dirpath + 'flux_i2.out'
     OPENR, 1, dirpath + 'out.tgyro.flux_i2'
     FOR ii = 0, N_it-1 DO BEGIN
        READF, 1, s
        READF, 1, s
        READF, 1, arr
        Gi2_neo[*,ii] = arr[1,*]
        Gi2_tur[*,ii] = arr[2,*]
        Qi2_neo[*,ii] = arr[3,*]
        Qi2_tur[*,ii] = arr[4,*]
     ENDFOR
     CLOSE, 1
  ENDIF

 ne_data_rho = -1
 ne_data = 0
 ne_data_err = 0
 OPENR, 1, dirpath + 'ne_data.txt', ERR=err
 IF (err EQ 0) THEN BEGIN
     READF, 1, n_data_pts
     arr = FLTARR(3,n_data_pts)
     READF, 1, arr
     CLOSE, 1

     ne_data_rho = REFORM(arr[0,*])
     ne_data = REFORM(arr[1,*])
     ne_data_err = REFORM(arr[2,*])
 ENDIF

 te_data_rho = -1
 te_data = 0
 te_data_err = 0
 OPENR, 1, dirpath + 'te_data.txt', ERR=err
 IF (err EQ 0) THEN BEGIN
     READF, 1, n_data_pts
     arr = FLTARR(3,n_data_pts)
     READF, 1, arr
     CLOSE, 1

     te_data_rho = REFORM(arr[0,*])
     te_data = REFORM(arr[1,*])
     te_data_err = REFORM(arr[2,*])
 ENDIF

 ti_data_rho = -1
 ti_data = 0
 ti_data_err = 0
 OPENR, 1, dirpath + 'ti_data.txt', ERR=err
 IF (err EQ 0) THEN BEGIN
     READF, 1, n_data_pts
     arr = FLTARR(3,n_data_pts)
     READF, 1, arr
     CLOSE, 1

     ti_data_rho = REFORM(arr[0,*])
     ti_data = REFORM(arr[1,*])
     ti_data_err = REFORM(arr[2,*])
 ENDIF

 data = {simdir: simdir, $     ;simulation directory

          ;experimental profile info, form FLTARR[EXP_N_RHO]
          exp_n_rho: exp_n_rho, $ ;# of datapoints in INPUT_profiles
          exp_rho: exp_rho, $   ;normalized toroidal flux
          exp_rmin: exp_rmin, $ ;r_min/a
          a: a, $               ;value of r_min on LCFS (m)
          exp_Rmaj: exp_Rmaj, $ ;Rmaj/a
          exp_q: exp_q, $       ;saftey factor
          exp_kappa: exp_kappa, $ ;elongation
          exp_delta: exp_delta, $ ;triangularity
          exp_zeta: exp_zeta, $ ;squareness
          exp_zmag: exp_zmag, $ ;Z0(r)  ;m
          exp_omega0: exp_omega0, $     ; 1/s
          exp_Te: exp_Te, $     ; keV
          exp_ne: exp_ne, $     ; 10**19/m**3
          exp_Zeff: exp_Zeff, $ 
          exp_flow_mom: exp_flow_mom, $ ;Nm
          exp_P_e: exp_Pe, $     ; total electron heating in Mw
          exp_P_i: exp_Pi, $     ; total ion heating in Mw
          exp_Pexch: exp_Pexch, $  ;ion-electron exchange term in Mw
          exp_flow_beam: exp_flow_beam, $ ; beam-driven Gamma_e  in MW/keV?
          exp_flow_wall: exp_flow_wall, $ ; wall-source-driven Gamma_e
          exp_flow: exp_flow, $  ; total exp. Gamma_e (beam+wall)
          exp_ni: exp_ni, $     ;main ion density, 10**19/m**3
          exp_ni2: exp_ni2, $   ;1st impurity density
	  exp_ni3: exp_ni3, $
          exp_Ti: exp_Ti, $     ;keV
          exp_Ti2: exp_Ti2, $
          exp_Ti3: exp_Ti3, $
          exp_Vtor: exp_Vtor, $  ;m/s
          exp_Vtor2: exp_Vtor2, $
          exp_Vtor3: exp_Vtor3, $
          exp_Vpol: exp_Vpol, $
          exp_Vpol2: exp_Vpol2, $
          exp_Vpol3: exp_Vpol3, $

          ;vector size for TGYRO output
          NX: NX, $             ;# of radial locations inc. 0
          N_it: N_it, $         ;total # of iterations

          ;following are arrays of form FLTARR[NX]
          rmin: rmin, $         ;rmin/a
          rho: rho, $           ;normalized toroidal flux
          q: q, $               ;flux-surface avg. saftey factor
          s_hat: s_hat, $       ;(r_min/q)*dq/dr_min
          kappa: kappa, $       ;plasma elongation
          s_kappa: s_kappa, $   ;(r_min/kappa)*dkappa/dr_min
          delta: delta, $       ;plasma triangularity
          s_delta: s_delta, $   ;(r_min/delta)*ddelta/dr_min
          shift: shift, $       ;Shafranov shift
          Rmaj: Rmaj, $         ;Rmaj/a
          b_unit: b_unit, $     ;b_unit (T)
	  z_mag: zmag, $        ;elevation of flux surfaces
	  dzmag: dzmag, $       ;d(z_mag)/dr
	  zeta: zeta, $         ;squareness
	  s_zeta: s_zeta, $     ;squareness shear
	  vol: vol, $           ;volume (m**3)
          Vprime: dvoldr, $     ;dV/dr_min, m**2

          ;following are arrays of form FLTARR[NX,N_it]
          Di_neo: Di_neo, $     ;Di_neo/chi_gB
          Di_tur: Di_tur, $     ;Di_turb/chi_gB
          De_neo: De_neo, $     ;De_neo/chi_gB
          De_tur: De_tur, $     ;De_turb/chi_gB
          chii_neo: chii_neo, $ ;chii_neo/chi_gB
          chii_tur: chii_tur, $ ;chii_turb/chi_gB
          chie_neo: chie_neo, $ ;chie_neo/chi_gB
          chie_tur: chie_tur, $ ;chie_turb/chi_gB
          Gi_neo: Gi_neo, $     ;particle fluxes in
          Gi_tur: Gi_tur, $     ;gB units
          Ge_neo: Ge_neo, $
          Ge_tur: Ge_tur, $
          Ge_tot: Ge_tot, $
          Ge_target: Ge_target, $
          Qi_neo: Qi_neo, $     ;energy fluxes in
          Qi_tur: Qi_tur, $     ;gB units
          Qi_tot: Qi_tot, $
          Qi_target: Qi_target, $
          Qe_neo: Qe_neo, $
          Qe_tur: Qe_tur, $
          Qe_tot: Qe_tot, $
          Qe_target: Qe_target, $
          Mflux_i_neo: Mflux_i_neo, $  ;momentum fluxes in
          Mflux_i_tur: Mflux_i_tur, $  ;gB units
          Mflux_e_neo: Mflux_i_neo, $
          Mflux_e_tur: Mflux_i_tur, $
          Mflux_tot: Mflux_tot, $
          Mflux_target: Mflux_target, $
          Expwd_i_tur: Expwd_i_tur, $ ;ion turbulent exchange
          Expwd_e_tur: Expwd_e_tur, $ ;elec turbulent exchange
          chi_gB: chi_gB, $      ;chi_gB in m**/s
          G_gB: G_gB, $          ;G_gB in 10**19/m**2/s
          Q_gB: Q_gB, $          ;Q_gB in MW/m**2
          Pi_gB: Pi_gB, $        ;Pi_gB in J/m**2
          c_s: c_s, $            ;c_s in m/s
          nu_ii: nu_ii, $        ;(a/c_s)/t_ii
          nu_ee: nu_ee, $        ;(a/c_s)/t_ee
          nu_exch: nu_exch, $    ;(a/c_s)nu_exch
          rho_i_star: rho_i_star, $ ;rho_i/a
          rho_star: rho_star, $  ;rho_s/a
          a_over_Lni: a_over_Lni, $
          a_over_Lne: a_over_Lne, $
          a_over_LTi: a_over_LTi, $
          a_over_LTe: a_over_LTe, $
          a_over_Lp: a_over_Lp, $
          gamma_e: gamma_e, $    ;a*gamma_e/c_s
          gamma_p: gamma_p, $    ;a*gamma_p/c_s
          n_i: n_i/1e13, $       ;converted to 10**19/m**3
          n_e: n_e/1e13, $
          Ti: Ti, $              ;keV
          Te: Te, $
          betae_unit: betae_unit, $
          M: M, $                ;Mach number R*omega/c_s
          P_alpha: p_alpha, $    ;all powers in MW
          P_brem: p_brem, $
          P_exch: p_exch, $
          P_i_aux: p_i_aux, $
          P_e_aux: p_e_aux, $
          P_i: p_i, $
          P_e: p_e, $

          N_ion: N_ion, $
          Di2_neo: Di2_neo, $
          Di2_tur: Di2_tur, $
          chii2_neo: chii2_neo, $
          chii2_tur: chii2_tur, $
          Gi2_neo: Gi2_neo, $
          Gi2_tur: Gi2_tur, $
          Qi2_neo: Qi2_neo, $
          Qi2_tur: Qi2_tur, $

          ;expt. data points
         ne_data_rho: ne_data_rho, $
         ne_data: ne_data, $
         ne_data_err: ne_data_err, $
         te_data_rho: te_data_rho, $
         te_data: te_data, $
         te_data_err: te_data_err, $
         ti_data_rho: ti_data_rho, $
         ti_data: ti_data, $
         ti_data_err: ti_data_err $
         }
  RETURN, data
END ;get_tgyro_loc_data
