FUNCTION make_synrefl_array, data, n_e, ITMIN = itmin, $
  ITMAX = itmax,  INTERP=interp, OMEGA0 = omega0, N_tor_frac = n_tor_frac, $
  Zloc = zloc, LR = lr, LZ=lz, DX=dx, GLOBAL_LOC = GLOBAL_loc, $
  SHOW_WINDOWS = show_windows
;
; C. Holland, UCSD
; v1.0 1-2-2009 Adapted from make_CECE08_synarrays.pro
; v1.1 8-5-2009 Updated to include a/cs, defaults set to 133363 2500 ms
; v2.0 7-26-2010: Updated to read units.out, Zloc, Lr, Lz, dx 
;                         n_tor_frac, global_loc keywords
;
; takes PostGYRO data  (for profile, time info) +
; n_e=[n_theta_plot,n_R,n_n,n_time] complex arrays, generates 
; 5 synthetic low-k reflectometry channel pairs at same radial
; locations.  Uses 2D gaussian spot function, similar to CECE, values
; from AEW 2009 work.
; 
; ITMIN/ITMAX specify starting/stopping time indexes
; INTERP does time interpolation needed b/c calc done in lab frame
; generally spacing = 0.25 should be ok.
; OMEGA0 is equilibrium ExB rotation frequency (get from run.out)
; GLOBAL_LOC = is array of positions for left-side channels w/ arb
; spacing, if not set, defaults to 5 equally spaced pairs
;  n_tor_Frac = # of equally spaced toroidal locations for 2D calculations
;
  DEFAULT, itmin, data.n_time/2
  DEFAULT, itmax, itmin
  sf = data.profile_data.theta_mult
  DEFAULT, interp, 1
  DEFAULT, omega0, 0
  DEFAULT, zloc, 5.5  ;cm
  DEFAULT, lr, 1.0 ;cm
  DEFAULT, lz, 3.5 ;cm
  DEFAULT, dx, 0.5 ;m
  DEFAULT, n_tor_frac, 1
  delta_tor_frac = 1./n_tor_frac

  OPENR, 1, GETENV('GYRO_DIR') + '/sim/' + data.simdir + '/units.out'
  READF, 1, tmp  ;m_ref
  READF, 1, tmp  ;b_unit
  READF, 1, Aphys
  Aphys *= 100 ; convert to cm
  PRINT, 'units.out: Aphys (cm) = ', Aphys
  READF, 1, tmp  ;c_s/a
  a_over_cs = 1e6/tmp  ;conver to microsec
  PRINT, 'units.out: a_over_cs (microsec) = ', a_over_cs
  CLOSE, 1

  ;create GYRO (R,Z) coords, set up triangles for interpolation
  ny = sf*data.n_theta_plot
  GYRO_GENERATE_RZCOORDS, data.profile_data, ny, GYRO_R, GYRO_Z

  ;create N_REFL *pairs* of REFL channels (for statistics) 
  ;spot = exp(-0.5*((R-R0)^2/Lr^2 + (Z-Z0)^2/Lz^2))
  ;when R-R0=2*Lr, or Z-Z0=2*Lz, spot = 1/e**2
  ;radial 1/e**2 diameter=4*L is 1 cm radially, 3.8 cm vertically

  IF N_ELEMENTS(global_loc) GT 0 THEN $
    N_REFL = N_ELEMENTS(global_loc) ELSE N_REFL = 5
  REFL_zloc = zloc/Aphys
  REFL_gyroidx = LONARR(2*N_REFL)
  REFL_spot = FLTARR(data.n_r,ny,2*N_REFL)
  REFL_spot_norm = FLTARR(2*N_REFL)
  REFL_lr = lr/(4*Aphys)
  REFL_lz = lz/(4.*Aphys)

  REFL_dx = dx/Aphys
  Lx = data.r[data.n_r-1]-data.r[0]
  IF N_ELEMENTS(global_loc) GT 0 THEN REFL_r0 = global_loc $
  ELSE REFL_r0 = data.r[0] + Lx*(FINDGEN(N_REFL)+1)/(N_REFL+1)
  REFL_xloc = FLTARR(2*N_REFL)
  REFL_rmid = FLTARR(2*N_REFL)
  for ir=0,N_REFL-1 do begin
      mind = MIN(ABS(REFL_r0[ir] - data.r),idx)
      REFL_rmid[2*ir] = data.profile_data.r[idx] 
      REFL_xloc[2*ir] = data.profile_data.R0[idx] + data.r[idx]
      mind = MIN(ABS(REFL_r0[ir]+REFL_dx - data.r),idx)
      REFL_rmid[2*ir+1] = data.profile_data.r[idx] 
      REFL_xloc[2*ir+1] = data.profile_data.R0[idx] + data.r[idx]
  ENDFOR

  FOR ir = 0, 2*N_REFL-1 DO BEGIN
      mind = MIN((GYRO_R - REFL_xloc[ir])^2 + (GYRO_Z - REFL_zloc)^2, $
                 minidx)
      REFL_gyroidx[ir] = minidx
      REFL_spot[*,*,ir] = EXP(-0.5*((GYRO_R - REFL_xloc[ir])^2/REFL_lr^2+$
                                    (GYRO_Z - REFL_zloc)^2/REFL_lz^2))
      REFL_spot_norm[ir] = TOTAL(REFL_spot[*,*,ir])
  ENDFOR


  IF KEYWORD_SET(show_windows) THEN BEGIN
      WINDOW, 0
      it = itmin
      ne_RZ = TRANSPOSE(GYRO_RTHETA_TRANSFORM(n_e[*,*,*,it], $
                                              data.profile_data,sf,$
                                              NU_SILENT=nu_silent, $
                                              TOR_FRAC=tor_frac))
      NR = data.n_r
      xr = [(data.profile_data.R0[0]+data.profile_data.r[0]),$
            (data.profile_data.R0[NR-1]+data.profile_data.r[NR-1])]*Aphys
      x0 =(data.profile_data.R0[NR/2]+data.profile_data.r[NR/2])*Aphys
      LX = 0.8*(xr[1] - xr[0])
      GYRO_RZ_COLOR_CONTOUR, TRANSPOSE(ne_RZ), GYRO_R, GYRO_Z, $
        XRANGE = [0.98*xr[0],1.02*xr[1]], $
        YRANGE=[Zloc-LX/2,Zloc+LX/2], $
        Aphys=Aphys, /XS, /YS, $
        TITLE = '!4d!Xn!De!N (t = ' + NUMTOSTRING(data.t[it]) + ')'
      OPLOT, GYRO_R[data.profile_data.n_bnd,*]*Aphys, $
             GYRO_Z[data.profile_data.n_bnd,*]*Aphys,$
             LINESTYLE=3
      OPLOT, GYRO_R[data.n_r-data.profile_data.n_bnd-1,*]*Aphys, $
             GYRO_Z[data.n_r-data.profile_data.n_bnd-1,*]*Aphys,$
             LINESTYLE=3

      FOR ir = 0, 2*N_REFL-1 DO BEGIN
          PLOTS, REFL_xloc[ir]*Aphys, REFL_zloc*Aphys, psym=4
          PLOTS, GYRO_R(REFL_gyroidx[ir])*Aphys, $
                 GYRO_Z(REFL_gyroidx[ir])*Aphys, $
                 psym=2, color=50
      ENDFOR
      CONTOUR, REFL_spot[*,*,N_REFL-1], GYRO_R*Aphys, GYRO_Z*Aphys, $
               LEVELS=[0.1,0.5,0.9], THICK=2,$
               C_LINESTYLE=[4,2,0], /OVERPLOT
      FOR ir=0, 2*N_REFL-1 DO $
             CONTOUR, REFL_spot[*,*,ir], GYRO_R*Aphys, GYRO_Z*Aphys, $
                      LEVELS=[0.1,0.5,0.9], THICK=2, $
                      C_LINESTYLE=[4,2,0], /OVERPLOT
     
  ENDIF

  ;set up time axis and result storage arrays
  NT0 = itmax - itmin + 1
  NT = 1 + interp*(NT0-1)
  t = itmin + FINDGEN(NT)/interp
  t_int = itmin+ INDGEN(NT)/interp
  eps = t - t_int

  GYRO_ne_sig = FLTARR(2*N_REFL, N_tor_frac, NT)
  REFL_sig = FLTARR(2*N_REFL,N_tor_frac, NT)
  C_I = COMPLEX(0,1)

  ;begin iterations
  t0=systime(1)
  FOR it = 0,NT-1 DO BEGIN
      elapsed=float(systime(1)-t0)
      eta=fix(elapsed/(it+1)*(NT-1-it))
      elapsed=fix(elapsed)
      hours=eta/3600
      hours_elapsed=elapsed/3600
      minutes=(eta-3600*hours)/60
      minutes_elapsed=(elapsed-3600*hours_elapsed)/60
      sec=(eta-3600*hours-60*minutes)
      seconds_elapsed=(elapsed-3600*hours_elapsed-60*minutes_elapsed)
      eta_str=string(hours,format='(I02)')+':'+string(minutes,format='(I02)')$
         +':'+string(sec,format='(I02)')
      elapsed_str=string(hours_elapsed,format='(I02)')+':'+$
         string(minutes_elapsed,format='(I02)')+':'+$
         string(seconds_elapsed,format='(I02)')
      PRINT, strtrim(it,2)+'/'+strtrim(NT-1,2),'     ETA:  '+eta_str+$
         '  Elapsed:  '+elapsed_str
         
      ;translate GYRO to rtheta space
      ne_lf = n_e[*,*,*,t_int[it]]

      ;do linear time interpolation
      IF (eps[it] GT 0) THEN BEGIN
          ne_lf2 = n_e[*,*,*,t_int[it]+1]
          ne_lf = (1. - eps[it])*ne_lf + eps[it]*ne_lf2
      ENDIF

      ;apply doppler shift
      FOR i_n = 1, data.n_n-1 DO BEGIN
          ne_lf[*,*,i_n] *= EXP(-C_I*data.n[i_n]*omega0*t[it])
      ENDFOR

      FOR i_tf = 0, n_tor_frac-1 DO BEGIN
          tor_frac = i_tf*delta_tor_frac
          ne_RZ = TRANSPOSE(GYRO_RTHETA_TRANSFORM(ne_lf, data.profile_data,sf,$
                                                  NU_SILENT=nu_silent, $
                                                  TOR_FRAC=tor_frac))
          
          nu_silent = 1
          FOR ir=0,2*N_REFL-1 DO BEGIN
              GYRO_ne_sig[ir,i_tf,it] = ne_RZ[REFL_gyroidx[ir]]
              REFL_sig[ir,i_tf,it] = TOTAL(ne_RZ*REFL_spot[*,*,ir])/$
                REFL_spot_norm[ir]
          ENDFOR
      ENDFOR
  ENDFOR


  IF KEYWORD_SET(show_windows) THEN BEGIN
      WINDOW, 1
      PLOT, t, GYRO_ne_sig[0,0,*], /XS, title='GYRO n!De!N, syn. REFL (red)'
      OPLOT, t, REFL_sig[0,0,*], COLOR=100
  ENDIF

  result = {SF:sf, NT:nt, T:t, a_over_cs: a_over_cs, INTERP:interp, $
            GYRO_ne:GYRO_ne_Sig, REFL:REFL_sig, $
            REFL_XLOC:REFL_xloc, REFL_ZLOC:REFL_zloc, $
            GYRO_REFL_XLOC:GYRO_R[REFL_gyroidx], $
            GYRO_REFL_ZLOC:GYRO_Z[REFL_gyroidx], $
            ITMIN:itmin, ITMAX:itmax, N_REFL:N_REFL, $
            REFL_rmid: REFL_rmid, $
            Aphys:Aphys, N_TOR_FRAC:N_tor_frac, DIR:data.simdir}

  RETURN, result
END ;make_synrefl_array
