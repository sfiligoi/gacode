FUNCTION read_input_profiles, FILENAME=filename
;
; C. Holland, UCSD
;
; v1.0: June 25 2010
;
; Reads in INPUT_profiles returns everything in one large IDL structure.
;
; FILENAME: full pathname; defaults to INPUT_profiles
;
; Usage: To load INPUT_profiles in a local directory, simply use
; IDL> data = get_tgyro_loc_data()
;

  IF N_ELEMENTS(filename) EQ 0 THEN filename = 'INPUT_profiles'
  
  ;first read INPUT_profiles to get inital info on fine grid
  s = STRING('#')
  OPENR, 1, filename
  WHILE((STRPOS(s,'#') EQ 0) OR (STRPOS(s,'V') EQ 0) OR (STRPOS(s,'O') EQ 0)) DO BEGIN
      READF, 1, s
  ENDWHILE
  n_rho = FIX(STRMID(s,6,3))
  READF, 1, s
  READF, 1, s
  arr = FLTARR(5,n_rho)

  WHILE (STRPOS(s, 'rho') EQ -1) DO READF, 1, s
  READF, 1, arr
  rho = REFORM(arr[0,*])
  rmin = REFORM(arr[1,*])
  a = rmin[n_rho-1]
  rmin /= a
  Rmaj = REFORM(arr[2,*])/a
  q = REFORM(arr[3,*])
  kappa = REFORM(arr[4,*])

  s = '#'
  WHILE (STRPOS(s, 'delta') EQ -1) DO READF, 1, s
  READF, 1, arr
  delta = REFORM(arr[0,*])
  Te = REFORM(arr[1,*])
  n_e = REFORM(arr[2,*])
  Zeff = REFORM(arr[3,*])
  Er = REFORM(arr[4,*])

  s = '#'
  WHILE (STRPOS(s, 'pow_e') EQ -1) DO READF, 1, s
  READF, 1, arr
  flow_mom = REFORM(arr[0,*])  ;Nm
  Pe = REFORM(arr[1,*])        ;Mw
  Pi = REFORM(arr[2,*])        ;Mw
  Pexch = REFORM(arr[3,*])     ;Mw
  zeta = REFORM(arr[4,*])

  READF, 1, s
  READF, 1, s
  READF, 1, arr
  flow_beam = REFORM(arr[0,*])  ;units?? MW/keV?
  flow_wall = REFORM(arr[1,*])
  flow = flow_beam + flow_wall
  zmag = REFORM(arr[3,*])/a ;Z elevation of flux surface, normalized to a
  ;columns 2-4 unused

  READF, 1, s
  READF, 1, s
  READF, 1, arr
  ni = REFORM(arr[0,*])    ;10**19/m**3
  ni2 = REFORM(arr[1,*])
  ni3 = REFORM(arr[2,*])
  ni4 = REFORM(arr[3,*])
  ni5 = REFORM(arr[4,*])

  READF, 1, s
  READF, 1, s
  READF, 1, arr
  Ti = REFORM(arr[0,*])
  Ti2 = REFORM(arr[1,*])
  Ti3 = REFORM(arr[2,*])
  Ti4 = REFORM(arr[3,*])
  Ti5 = REFORM(arr[4,*])

  READF, 1, s
  READF, 1, s
  READF, 1, arr
  Vtor = REFORM(arr[0,*])
  Vtor2 = REFORM(arr[0,*])
  Vtor3 = REFORM(arr[0,*])
  Vtor4 = REFORM(arr[0,*])
  Vtor5 = REFORM(arr[0,*])

  READF, 1, s
  READF, 1, s
  READF, 1, arr
  Vpol = REFORM(arr[0,*])
  Vpol2 = REFORM(arr[0,*])
  Vpol3 = REFORM(arr[0,*])
  Vpol4 = REFORM(arr[0,*])
  Vpol5 = REFORM(arr[0,*])

  CLOSE, 1

 data = {filname: filename, $
          n_rho: n_rho, $ ;# of datapoints in INPUT_profiles
          rho: rho, $   ;normalized toroidal flux
          rmin: rmin, $ ;r_min/a
          a: a, $               ;value of r_min on LCFS (m)
          Rmaj: Rmaj, $ ;Rmaj/a
          q: q, $       ;saftey factor
          kappa: kappa, $ ;elongation
          delta: delta, $ ;triangularity
          zeta: zeta, $ ;squareness
          zmag: zmag, $ ;elevation (normed to a)
          Er: Er, $     ; kV/m  
          Te: Te, $     ; keV
          n_e: n_e, $     ; 10**19/m**3
          Zeff: Zeff, $ 
          flow_mom: flow_mom, $ ;Nm
          P_e: Pe, $     ; total electron heating in Mw
          P_i: Pi, $     ; total ion heating in Mw
          Pexch: Pexch, $  ;ion-electron exchange term in Mw
          flow_beam: flow_beam, $ ; beam-driven Gamma_e  in MW/keV?
          flow_wall: flow_wall, $ ; wall-source-driven Gamma_e
          flow: flow, $  ; total exp. Gamma_e (beam+wall)
          ni: ni, $     ;main ion density, 10**19/m**3
          ni2: ni2, $   ;1st impurity density
          ni3: ni3, $
          Ti: Ti, $     ;keV
          Ti2: Ti2, $
          Ti3: Ti3, $
          Vtor: Vtor, $  ;m/s
          Vtor2: Vtor2, $
          Vtor3: Vtor3, $
          Vpol: Vpol, $
          Vpol2: Vpol2, $
          Vpol3: Vpol3 }
  RETURN, data
END ;read_input_profiles
