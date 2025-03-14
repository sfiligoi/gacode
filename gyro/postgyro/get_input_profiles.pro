FUNCTION get_input_profiles, simdir, FILENAME=filename, $
  EXTRA_FILENAME=extra_filename, DIRLOC=dirloc
;
; C. Holland, UCSD
;
; v1.0: Sept. 15,2011 stolen from get_loc_tgyro_data.pro
; v2.0: May 11, 2012 updated to read in input.profiles files with arbitrary names
;	(e.g. input.profiles.orig, old.input.prof, etc.)
; v2.1: March 13, 2013 updated to get vprime for input.profiles.extra
; v2.2: July 19, 2013 added Bt_exp and arho_exp to output structure
; v3.0: May 24, 2016 updated for current GACODE input.profiles structure
;
; Reads in input.profiles from a local subdirectory
;
; KEYWORDS
; simdir: string containing name of valid subdirectory in local directory
; filename: name of input.profiles file to open, defaults to
; input.profiles
; extra_filename: name of input.profiles.extra filename to open,
; defaults to input.profiles.extra
; dirloc: string with pathname to directory containing simdir
;
  IF N_ELEMENTS(simdir) EQ 0 THEN BEGIN
      MESSAGE, 'Need to specify a directory!', /INFO
      RETURN, 0
  ENDIF

  DEFAULT, dirloc, '.'
  dirpath = dirloc + '/' + simdir + '/'

  DEFAULT, filename, 'input.profiles'
  DEFAULT, extra_filename, 'input.profiles.extra'

  s = STRING('#')
  OPENR, 1, dirpath + filename, ERROR=err
  IF (err NE 0) THEN BEGIN
      PRINT, "Couldn't open " + dirpath + ' ' + filename + "!  Returning 0"
      RETURN, 0
  ENDIF

  WHILE((STRPOS(s,'#') EQ 0) OR (STRPOS(s,'V') EQ 0) OR (STRPOS(s,'O') EQ 0)) DO BEGIN
      READF, 1, s
  ENDWHILE
  exp_n_ion = FIX(STRMID(s,6))
  READF, 1, s
  exp_n_rho = FIX(STRMID(s,6))
  READF, 1, s
  BT_exp = DOUBLE(STRMID(s,7))
  READF, 1, s
  arho_exp = DOUBLE(STRMID(s,9))

  arr = DBLARR(5,exp_n_rho)
  WHILE (STRPOS(s, 'rho') EQ -1) DO READF, 1, s
  READF, 1, arr
  exp_rho = REFORM(arr[0,*])
  exp_rmin = REFORM(arr[1,*])
  a = exp_rmin[exp_n_rho-1]
  exp_rmin /= a
  exp_psi = REFORM(arr[2,*])
  exp_psiN = (exp_psi - exp_psi[0])/(exp_psi[exp_n_rho-1] - exp_psi[0])
  exp_q = REFORM(arr[3,*])
  exp_omega0 = REFORM(arr[4,*])

  s = '#'
  WHILE (STRPOS(s, 'rmaj') EQ -1) DO READF, 1, s
  READF, 1, arr
  exp_rmaj = REFORM(arr[0,*])/a
  exp_zmag = REFORM(arr[1,*])
  exp_kappa = REFORM(arr[2,*])
  exp_delta = REFORM(arr[3,*])
  exp_zeta = REFORM(arr[4,*])

  s = '#'
  WHILE (STRPOS(s, 'ne') EQ -1) DO READF, 1, s
  READF, 1, arr
  exp_ne = REFORM(arr[0,*])/a
  exp_Te = REFORM(arr[1,*])
  exp_ptot = REFORM(arr[2,*])
  exp_zeff = REFORM(arr[3,*])
  null = REFORM(arr[4,*])

  s = '#'
  WHILE (STRPOS(s, 'ni_1') EQ -1) DO READF, 1, s
  READF, 1, arr
  exp_ni1 = REFORM(arr[0,*])
  exp_ni2 = REFORM(arr[1,*])
  exp_ni3 = REFORM(arr[2,*])
  exp_ni4 = REFORM(arr[3,*])
  exp_ni5 = REFORM(arr[4,*])
  s = '#'
  WHILE (STRPOS(s, 'ni_6') EQ -1) DO READF, 1, s
  READF, 1, arr
  exp_ni6 = REFORM(arr[0,*])
  exp_ni7 = REFORM(arr[1,*])
  exp_ni8 = REFORM(arr[2,*])
  exp_ni9 = REFORM(arr[3,*])
  exp_ni10 = REFORM(arr[4,*])

  s = '#'
  WHILE (STRPOS(s, 'Ti_1') EQ -1) DO READF, 1, s
  READF, 1, arr
  exp_Ti1 = REFORM(arr[0,*])
  exp_Ti2 = REFORM(arr[1,*])
  exp_Ti3 = REFORM(arr[2,*])
  exp_Ti4 = REFORM(arr[3,*])
  exp_Ti5 = REFORM(arr[4,*])
  s = '#'
  WHILE (STRPOS(s, 'Ti_6') EQ -1) DO READF, 1, s
  READF, 1, arr
  exp_Ti6 = REFORM(arr[0,*])
  exp_Ti7 = REFORM(arr[1,*])
  exp_Ti8 = REFORM(arr[2,*])
  exp_Ti9 = REFORM(arr[3,*])
  exp_Ti10 = REFORM(arr[4,*])

  s = '#'
  WHILE (STRPOS(s, 'vtor_1') EQ -1) DO READF, 1, s
  READF, 1, arr
  exp_vtor1 = REFORM(arr[0,*])
  exp_vtor2 = REFORM(arr[1,*])
  exp_vtor3 = REFORM(arr[2,*])
  exp_vtor4 = REFORM(arr[3,*])
  exp_vtor5 = REFORM(arr[4,*])
  s = '#'
  WHILE (STRPOS(s, 'vtor_6') EQ -1) DO READF, 1, s
  READF, 1, arr
  exp_vtor6 = REFORM(arr[0,*])
  exp_vtor7 = REFORM(arr[1,*])
  exp_vtor8 = REFORM(arr[2,*])
  exp_vtor9 = REFORM(arr[3,*])
  exp_vtor10 = REFORM(arr[4,*])

  s = '#'
  WHILE (STRPOS(s, 'vpol_1') EQ -1) DO READF, 1, s
  READF, 1, arr
  exp_vpol1 = REFORM(arr[0,*])
  exp_vpol2 = REFORM(arr[1,*])
  exp_vpol3 = REFORM(arr[2,*])
  exp_vpol4 = REFORM(arr[3,*])
  exp_vpol5 = REFORM(arr[4,*])
  s = '#'
  WHILE (STRPOS(s, 'vpol_6') EQ -1) DO READF, 1, s
  READF, 1, arr
  exp_vpol6 = REFORM(arr[0,*])
  exp_vpol7 = REFORM(arr[1,*])
  exp_vpol8 = REFORM(arr[2,*])
  exp_vpol9 = REFORM(arr[3,*])
  exp_vpol10 = REFORM(arr[4,*])

  s = '#'
  WHILE (STRPOS(s, 'flow_beam') EQ -1) DO READF, 1, s
  READF, 1, arr
  exp_flow_beam = REFORM(arr[0,*])  ;kW/eV
  exp_flow_wall = REFORM(arr[1,*])        ;kW/ev
  exp_flow_mom = REFORM(arr[2,*])        ;Nm
  null = REFORM(arr[3,*])
  null = REFORM(arr[4,*])

  s = '#'
  WHILE (STRPOS(s, 'pow_e') EQ -1) DO READF, 1, s
  READF, 1, arr
  exp_pow_e = REFORM(arr[0,*])  ;MW
  exp_pow_i = REFORM(arr[1,*])        ;MW
  exp_pow_exch = REFORM(arr[2,*])        ;MW
  exp_pow_e_aux = REFORM(arr[3,*])   ;MW
  exp_pow_i_aux = REFORM(arr[4,*])   ;MW

  s = '#'
  WHILE (STRPOS(s, 'pow_e_fus') EQ -1) DO READF, 1, s
  READF, 1, arr
  exp_pow_e_fus = REFORM(arr[0,*])  ;MW
  exp_pow_i_fus = REFORM(arr[1,*])        ;MW
  exp_pow_e_sync = REFORM(arr[2,*])        ;MW
  exp_pow_e_brem = REFORM(arr[3,*])   ;MW
  exp_pow_e_line = REFORM(arr[4,*])   ;MW

  CLOSE, 1

  n_exp_profile = 46           ;should be 37, but use 25 for back-compatibility
  arr = FLTARR(exp_n_rho, n_exp_profile)
  OPENR, 1, dirpath + extra_filename, ERR=i_err
  IF (i_err NE 0) THEN PRINT, "Could not find " + extra_filename ELSE BEGIN
      ;need to skip over comment lines
      s = ' '
      n_cmt = 0
      readf, 1, s
      while (strpos(s,'#') EQ 0) do begin
          n_cmt += 1
          readf, 1, s
      endwhile
      close, 1
      openr, 1, dirpath + extra_filename
      for ii = 0, n_cmt-1 do readf, 1, s
      
      READF, 1, arr
      CLOSE, 1
  ENDELSE
  exp_Bunit = REFORM(arr[*,0])
  exp_vol = REFORM(arr[*,32])
  exp_vprime = REFORM(arr[*,33])

 data = {simdir: simdir, $     ;simulation directory

          ;experimental profile info, form FLTARR[EXP_N_RHO]
          exp_n_ion: exp_n_ion, $ ;# number of ions
          exp_n_rho: exp_n_rho, $ ;# of datapoints in INPUT_profiles
          BT_exp: BT_exp, $  ;on-axis Btor in T
          arho_exp: arho_exp, $ ;value of rho at rmin/a = 1 in m
          exp_rho: exp_rho, $   ;normalized toroidal flux
          exp_rmin: exp_rmin, $ ;r_min/a
          exp_polflux: exp_psi, $ ;poloidal flux (Wb/rad)
          exp_psiN: exp_psiN, $ ;normalized pol. flux
          a: a, $               ;value of r_min on LCFS (m)
          exp_rmaj: exp_rmaj, $ ;Rmaj/a
          exp_q: exp_q, $       ;saftey factor
          exp_kappa: exp_kappa, $ ;elongation
          exp_delta: exp_delta, $ ;triangularity
          exp_zeta: exp_zeta, $ ;squareness
          exp_zmag: exp_zmag, $ ;Z0(r)  ;m
          exp_vol: exp_vol, $ ;m^3
          exp_vprime: exp_vprime, $ ;dV/drmin m^2
          exp_Bunit: exp_Bunit, $ ;T
          exp_omega0: exp_omega0, $     ; 1/s
          exp_Te: exp_Te, $     ; keV
          exp_ne: exp_ne, $     ; 10**19/m**3
          exp_Zeff: exp_Zeff, $ 
          exp_flow_mom: exp_flow_mom, $ ;Nm
          exp_pow_e: exp_pow_e, $     ; total electron heating in Mw
          exp_pow_i: exp_pow_i, $     ; total ion heating in Mw
          exp_pow_exch: exp_pow_exch, $  ;ion-electron exchange term in Mw
          exp_pow_e_aux: exp_pow_e_aux, $     ; total electron heating in Mw
          exp_pow_i_aux: exp_pow_i_aux, $     ; total ion heating in Mw
  	  exp_pow_e_fus: exp_pow_e_fus, $    ; integrated e- heating by alphas
  	  exp_pow_i_fus: exp_pow_i_fus, $    ; integrated ionheating by alphas
  	  exp_pow_e_sync: exp_pow_e_sync, $  ; integrated sync. radiation
  	  exp_pow_e_brem: exp_pow_e_brem, $  ; integrated bremstrahlung radiation
  	  exp_pow_e_line: exp_pow_e_line, $  ; integrated line radiation
          exp_flow_beam: exp_flow_beam, $ ; beam-driven Gamma_e  in MW/keV
          exp_flow_wall: exp_flow_wall, $ ; wall-source-driven Gamma_e
          exp_flow: exp_flow_Beam+exp_flow_wall, $  ; total exp. Gamma_e (beam+wall)
          exp_ptot: exp_ptot, $ ;total pressure inc. fast ions (pa)
          exp_ni: exp_ni1, $     ;main ion density, 10**19/m**3
          exp_ni2: exp_ni2, $   ;1st impurity density
          exp_ni3: exp_ni3, $
          exp_ni4: exp_ni4, $
          exp_ni5: exp_ni5, $
          exp_ni6: exp_ni6, $
          exp_ni7: exp_ni7, $
          exp_ni8: exp_ni8, $
          exp_ni9: exp_ni9, $
          exp_ni10: exp_ni10, $
          exp_Ti: exp_Ti1, $     ;keV
          exp_Ti2: exp_Ti2, $
          exp_Ti3: exp_Ti3, $
          exp_Ti4: exp_Ti4, $
          exp_Ti5: exp_Ti5, $
          exp_Ti6: exp_Ti6, $
          exp_Ti7: exp_Ti7, $
          exp_Ti8: exp_Ti8, $
          exp_Ti9: exp_Ti9, $
          exp_Ti10: exp_Ti10, $
          exp_Vtor: exp_Vtor1, $  ;m/s
          exp_Vtor2: exp_Vtor2, $
          exp_Vtor3: exp_Vtor3, $
          exp_Vtor4: exp_Vtor4, $
          exp_Vtor5: exp_Vtor5, $
          exp_Vtor6: exp_Vtor6, $
          exp_Vtor7: exp_Vtor7, $
          exp_Vtor8: exp_Vtor8, $
          exp_Vtor9: exp_Vtor9, $
          exp_Vtor10: exp_Vtor10, $
          exp_Vpol: exp_Vpol1, $
          exp_Vpol2: exp_Vpol2, $
          exp_Vpol3: exp_Vpol3, $
          exp_Vpol4: exp_Vpol4, $
          exp_Vpol5: exp_Vpol5, $
          exp_Vpol6: exp_Vpol6, $
          exp_Vpol7: exp_Vpol7, $
          exp_Vpol8: exp_Vpol8, $
          exp_Vpol9: exp_Vpol9, $
          exp_Vpol10: exp_Vpol10 }

  RETURN, data
END ;get_input_profiles
