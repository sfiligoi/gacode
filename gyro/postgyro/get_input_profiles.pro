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
  exp_n_rho = FIX(STRMID(s,6,3))
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
  exp_polflux = REFORM(arr[4,*]) ;Wb/rad

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

  n_exp_profile = 25            ;should be 37, but use 25 for back-compatibility
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
  exp_vprime = REFORM(arr[*,23])

 data = {simdir: simdir, $     ;simulation directory

          ;experimental profile info, form FLTARR[EXP_N_RHO]
          exp_n_rho: exp_n_rho, $ ;# of datapoints in INPUT_profiles
          BT_exp: BT_exp, $  ;on-axis Btor in T
          arho_exp: arho_exp, $ ;value of rho at rmin/a = 1 in m
          exp_rho: exp_rho, $   ;normalized toroidal flux
          exp_rmin: exp_rmin, $ ;r_min/a
          a: a, $               ;value of r_min on LCFS (m)
          exp_Rmaj: exp_Rmaj, $ ;Rmaj/a
          exp_q: exp_q, $       ;saftey factor
          exp_kappa: exp_kappa, $ ;elongation
          exp_delta: exp_delta, $ ;triangularity
          exp_zeta: exp_zeta, $ ;squareness
          exp_zmag: exp_zmag, $ ;Z0(r)  ;m
          exp_vprime: exp_vprime, $ ;dV/drmin m^2
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
          exp_ptot: exp_ptot, $ ;total pressure inc. fast ions (pa)
          exp_polflux: exp_polflux, $ ;poloidal flux (Wb/rad)
          exp_ni: exp_ni, $     ;main ion density, 10**19/m**3
          exp_ni2: exp_ni2, $   ;1st impurity density
	  exp_ni3: exp_ni3, $
          exp_Ti: exp_Ti, $     ;keV
          exp_Ti2: exp_Ti2, $
          exp_Ti3: exp_Ti3, $
          exp_Ti4: exp_Ti4, $
          exp_Ti5: exp_Ti5, $
          exp_Vtor: exp_Vtor, $  ;m/s
          exp_Vtor2: exp_Vtor2, $
          exp_Vtor3: exp_Vtor3, $
          exp_Vtor4: exp_Vtor4, $
          exp_Vtor5: exp_Vtor5, $
          exp_Vpol: exp_Vpol, $
          exp_Vpol2: exp_Vpol2, $
          exp_Vpol3: exp_Vpol3, $
          exp_Vpol4: exp_Vpol4, $
          exp_Vpol5: exp_Vpol5 }

  RETURN, data
END ;get_input_profiles
