FUNCTION make_synCECE_array, data, Te, SF = sf, ITMIN = itmin, $
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
; Te=[n_theta_plot,n_R,n_n,n_time] complex arrays, generates 
; 5 synthetic CECE channel pairs at same radial locations.
; 
; SF is intepolation in theta (set equal to
; data.profile_data.n_theta_mult to use GYRO-caluclated nu_geo)
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
  DEFAULT, sf, data.profile_data.theta_mult
  DEFAULT, interp, 1
  DEFAULT, omega0, 0
  DEFAULT, zloc, 5.5  ;cm
  DEFAULT, lr, 1.1 ;cm
  DEFAULT, lz, 3.2 ;cm
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

  ;create N_CECE *pairs* of CECE channels (for statistics) at
  ;same locations as BES, but above rather than below midplane
  ;spot = exp(-0.5*((R-R0)^2/Lr^2 + (Z-Z0)^2/Lz^2))
  ;when R-R0=2*Lr, or Z-Z0=2*Lz, spot = 1/e**2
  ;radial 1/e**2 diameter=4*L is 1 cm radially, 3.8 cm vertically

  IF N_ELEMENTS(global_loc) GT 0 THEN $
    N_CECE = N_ELEMENTS(global_loc) ELSE N_CECE = 5
  CECE_zloc = zloc/Aphys
  CECE_gyroidx = LONARR(2*N_CECE)
  CECE_spot = FLTARR(data.n_r,ny,2*N_CECE)
  CECE_spot_norm = FLTARR(2*N_CECE)
  CECE_lr = lr/(4*Aphys)
  CECE_lz = lz/(4.*Aphys)

  CECE_dx = dx/Aphys
  Lx = data.r[data.n_r-1]-data.r[0]
  IF N_ELEMENTS(global_loc) GT 0 THEN CECE_r0 = global_loc $
  ELSE CECE_r0 = data.r[0] + Lx*(FINDGEN(N_CECE)+1)/(N_CECE+1)
  CECE_xloc = FLTARR(2*N_CECE)
  CECE_rmid = FLTARR(2*N_CECE)
  Te0 = FLTARR(2*N_CECE)
  for ir=0,N_CECE-1 do begin
      mind = MIN(ABS(CECE_r0[ir] - data.r),idx)
      CECE_rmid[2*ir] = data.profile_data.r[idx] 
      CECE_xloc[2*ir] = data.profile_data.R0[idx] + data.r[idx]
      Te0[2*ir] = data.profile_data.T[data.profile_data.n_spec-1,idx]
      mind = MIN(ABS(CECE_r0[ir]+CECE_dx - data.r),idx)
      CECE_rmid[2*ir+1] = data.profile_data.r[idx] 
      CECE_xloc[2*ir+1] = data.profile_data.R0[idx] + data.r[idx]
      Te0[2*ir+1] = data.profile_data.T[data.profile_data.n_spec-1,idx]
  ENDFOR

  FOR ir = 0, 2*N_CECE-1 DO BEGIN
      mind = MIN((GYRO_R - CECE_xloc[ir])^2 + (GYRO_Z - CECE_zloc)^2, $
                 minidx)
      CECE_gyroidx[ir] = minidx
      CECE_spot[*,*,ir] = EXP(-0.5*((GYRO_R - CECE_xloc[ir])^2/CECE_lr^2+$
                                    (GYRO_Z - CECE_zloc)^2/CECE_lz^2))
      CECE_spot_norm[ir] = TOTAL(CECE_spot[*,*,ir])
  ENDFOR


  IF KEYWORD_SET(show_windows) THEN BEGIN
      WINDOW, 0
      it = itmin
      Te_RZ = TRANSPOSE(GYRO_RTHETA_TRANSFORM(Te[*,*,*,it], $
                                              data.profile_data,sf,$
                                              NU_SILENT=nu_silent, $
                                              TOR_FRAC=tor_frac))
      NR = data.n_r
      xr = [(data.profile_data.R0[0]+data.profile_data.r[0]),$
            (data.profile_data.R0[NR-1]+data.profile_data.r[NR-1])]*Aphys
      x0 =(data.profile_data.R0[NR/2]+data.profile_data.r[NR/2])*Aphys
      LX = 0.8*(xr[1] - xr[0])
      GYRO_RZ_COLOR_CONTOUR, TRANSPOSE(te_RZ), GYRO_R, GYRO_Z, $
        XRANGE = [0.98*xr[0],1.02*xr[1]], $
        YRANGE=[Zloc-LX/2,Zloc+LX/2], $
        Aphys=Aphys, /XS, /YS, $
        TITLE = '!4d!XT!De!N (t = ' + NUMTOSTRING(data.t[it]) + ')'
      OPLOT, GYRO_R[data.profile_data.n_bnd,*]*Aphys, $
             GYRO_Z[data.profile_data.n_bnd,*]*Aphys,$
             LINESTYLE=3
      OPLOT, GYRO_R[data.n_r-data.profile_data.n_bnd-1,*]*Aphys, $
             GYRO_Z[data.n_r-data.profile_data.n_bnd-1,*]*Aphys,$
             LINESTYLE=3

      FOR ir = 0, 2*N_CECE-1 DO BEGIN
          PLOTS, CECE_xloc[ir]*Aphys, CECE_zloc*Aphys, psym=4
          PLOTS, GYRO_R(CECE_gyroidx[ir])*Aphys, $
                 GYRO_Z(CECE_gyroidx[ir])*Aphys, $
                 psym=2, color=50
      ENDFOR
      CONTOUR, CECE_spot[*,*,N_CECE-1], GYRO_R*Aphys, GYRO_Z*Aphys, $
               LEVELS=[0.1,0.5,0.9], THICK=2,$
               C_LINESTYLE=[4,2,0], /OVERPLOT
      FOR ir=0, 2*N_CECE-1 DO $
             CONTOUR, CECE_spot[*,*,ir], GYRO_R*Aphys, GYRO_Z*Aphys, $
                      LEVELS=[0.1,0.5,0.9], THICK=2, $
                      C_LINESTYLE=[4,2,0], /OVERPLOT
     
  ENDIF

  ;set up time axis and result storage arrays
  NT0 = itmax - itmin + 1
  NT = 1 + interp*(NT0-1)
  t = itmin + FINDGEN(NT)/interp
  t_int = itmin+ INDGEN(NT)/interp
  eps = t - t_int

  GYRO_Te_sig = FLTARR(2*N_CECE, N_tor_frac, NT)
  CECE_sig = FLTARR(2*N_CECE,N_tor_frac, NT)
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
      Te_lf = te[*,*,*,t_int[it]]

      ;do linear time interpolation
      IF (eps[it] GT 0) THEN BEGIN
          Te_lf2 = te[*,*,*,t_int[it]+1]
          Te_lf = (1. - eps[it])*Te_lf + eps[it]*Te_lf2
      ENDIF

      ;apply doppler shift
      FOR i_n = 1, data.n_n-1 DO BEGIN
          Te_lf[*,*,i_n] *= EXP(-C_I*data.n[i_n]*omega0*t[it])
      ENDFOR

      FOR i_tf = 0, n_tor_frac-1 DO BEGIN
          tor_frac = i_tf*delta_tor_frac
          Te_RZ = TRANSPOSE(GYRO_RTHETA_TRANSFORM(te_lf, data.profile_data,sf,$
                                                  NU_SILENT=nu_silent, $
                                                  TOR_FRAC=tor_frac))
          
          nu_silent = 1
          FOR ir=0,2*N_CECE-1 DO BEGIN
              GYRO_Te_sig[ir,i_tf,it] = Te_RZ[CECE_gyroidx[ir]]
              CECE_sig[ir,i_tf,it] = TOTAL(Te_RZ*CECE_spot[*,*,ir])/$
                CECE_spot_norm[ir]
          ENDFOR
      ENDFOR
  ENDFOR


  IF KEYWORD_SET(show_windows) THEN BEGIN
      WINDOW, 1
      PLOT, t, GYRO_Te_sig[0,0,*], /XS, title='GYRO Te, syn. CECE (red)'
      OPLOT, t, CECE_sig[0,0,*], COLOR=100
  ENDIF

  result = {SF:sf, NT:nt, T:t, a_over_cs: a_over_cs, INTERP:interp, $
            GYRO_Te:GYRO_Te_Sig, CECE:CECE_sig, $
            CECE_XLOC:CECE_xloc, CECE_ZLOC:CECE_zloc, $
            GYRO_CECE_XLOC:GYRO_R[CECE_gyroidx], $
            GYRO_CECE_ZLOC:GYRO_Z[CECE_gyroidx], $
            ITMIN:itmin, ITMAX:itmax, N_CECE:N_CECE, $
            CECE_rmid: CECE_rmid, Te0:Te0, $
            Aphys:Aphys, N_TOR_FRAC:N_tor_frac, DIR:data.simdir}

  RETURN, result
END ;make_syncece_array
