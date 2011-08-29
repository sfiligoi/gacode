FUNCTION make_synbes_array, data, n_e, TOR_FRAC=tor_frac, $
  ITMIN = itmin, ITMAX = itmax,  INTERP=interp, OMEGA0=omega0, $
  N_TOR_FRAC=n_tor_frac, NR_BESchannels = NR_BESchannels, $
  NZ_BESchannels = NZ_BESchannels, BES_Z00 = BES_Z00, BES_DR=BES_DR, $
  BES_DZ=BES_DZ, SHOW_WINDOWS=show_windows, $
  Aphys = Aphys, A_over_Cs = a_over_cs
;
; C. Holland, UCSD
; v1.0: 9.14.2010
; v1.1: 10.12.2010: added Aphys, a_over_cs keywords to enable external
; inputs for these parameters when using tgyro2gyro input files
; v2.0: 8.25.2011: updated for consistency with gacode file structures
;
; takes PostGYRO data  (for profile, time info) +
; ne=[n_theta_plot,n_r,n_n,n_time] complex arrays, generates 
; NRxNZxn_tor_fac synthetic BES channels
; 
; SF is intepolation in theta (set equal to
; data.profile_data.n_theta_mult to use GYRO-caluclated nu_geo)
; ITMIN/ITMAX specify starting/stopping time indexes
; INTERP does time interpolation needed b/c calc done in lab frame
; generally spacing = 0.25 should be ok.
; OMEGA0 is equilibrium ExB rotation frequency at mid-radius norm point)
;  n_tor_Frac = # of equally spaced toroidal locations for 2D
;  calculations
;  BES parameters specifcy # of (R,Z) points, and spacing between them
;

  IF KEYWORD_SET(show_windows) THEN windows_flag=1 ELSE windows_flag=0
  PRINT, 'windows_flag = ', windows_flag

  DEFAULT, omega0, data.w0[data.n_r/2]
  PRINT, 'omega0 (c_s/a): ', omega0
  DEFAULT, n_tor_frac, 1
  PRINT, '# toroidal locations: ', n_tor_frac
  delta_tor_frac = 1./n_tor_frac;
  sf = data.theta_mult
  PRINT, 'sf = ', sf
  DEFAULT, itmin, data.n_time/2
  PRINT, 'itmin = ', itmin
  DEFAULT, itmax, data.n_time-1
  PRINT, 'itmax = ', itmax
  DEFAULT, interp, 1
  PRINT, 'interp = ', interp
  DEFAULT, NR_BESchannels, 5
  PRINT, 'NR_BESchannels = ', NR_BESchannels
  DEFAULT, NZ_BESchannels, 6
  PRINT, 'NZ_BESchannels = ', NZ_BESchannels
  DEFAULT, BES_Z00, -4. ;cm
  PRINT, 'BES Z00 (cm): ', BES_Z00
  DEFAULT, BES_DR, 0.9 ;cm
  PRINT, 'BES DR (cm): ', BES_DR
  DEFAULT, BES_DZ, 1.2 ;cm
  PRINT, 'BES DZ (cm): ', BES_DZ

  DEFAULT, Aphys, data.Aphys
  PRINT, 'Aphys (cm) = ', Aphys
  DEFAULT, a_over_cs, 1.e3/data.csda  ;data.csda in kHz
  PRINT, 'a_over_cs (microsec) = ', a_over_cs

  ;create GYRO (R,Z) coords, set up triangles for interpolation
  ny = sf*data.n_theta_plot
  GYRO_GENERATE_RZCOORDS, data, ny, GYRO_R, GYRO_Z

  ;start with BES
;  psf = psf_in   ;will alter some fields in psf, so work on temp
;  compy
  RESTORE, GETENV('GYRO_DIR') + '/sim/' + data.simdir + '/psf.sav'
  psf.r = (psf.r - psf.stats.r)/Aphys
  psf.p = (psf.p - psf.stats.z)/Aphys

  psf_NR= N_ELEMENTS(psf.r)
  psf_NZ = N_ELEMENTS(psf.p)
  psf_dx = (psf.r[psf_NR-1] - psf.r[0])/(psf_NR-1)
  psf_dz = (psf.p[psf_NZ-1] - psf.p[0])/(psf_NZ-1)

  psf_interp = FLTARR(data.n_r, ny, NR_BESchannels, NZ_BESchannels)
  psf_interp_norm = FLTARR(NR_BESchannels, NZ_BESchannels)
  BES_xloc = FLTARR(NR_BESchannels, NZ_BESchannels)
  BES_zloc = FLTARR(NR_BESchannels, NZ_BESchannels)
  BES_gyroidx = LONARR(NR_BESchannels, NZ_BESchannels)

  ;R00 = major radius/a of bes channel 9 is at middle of shot
  ir = data.n_r/2
  R00 =  (data.R0[ir] + data.r[ir])
  Z00 = BES_Z00/Aphys

  ;create array BES psfs interpolated onto GYRO grid
  FOR ir = 0, NR_BESchannels-1 DO FOR iz = 0, NZ_BESchannels-1 DO BEGIN
      BES_xloc[ir,iz] = R00 + (ir-NR_BESchannels/2)*BES_DR/Aphys
      BES_zloc[ir,iz] = Z00 - iz*BES_DZ/Aphys

      ;for now, use minidx to find closest GYRO point; maybe interpolate later
      mind = MIN((GYRO_R - BES_xloc[ir,iz])^2 + (GYRO_Z - BES_zloc[ir,iz])^2, $
                 minidx)
      BES_gyroidx[ir,iz] = minidx

      ;regrid PSF onto GYRO coords
      ir_GYRO = (GYRO_R - psf.R[0] - BES_xloc[ir,iz])/psf_dx
      iz_GYRO = (GYRO_Z - psf.P[0] - BES_zloc[ir,iz])/psf_dz
      psf_interp[*,*,ir,iz] = BILINEAR(psf.psf, ir_GYRO, iz_GYRO, MISSING=0)
      psf_interp_norm[ir,iz] = TOTAL(psf_interp[*,*,ir,iz])
  ENDFOR

  it = itmin
  ne_RZ = TRANSPOSE(GYRO_RTHETA_TRANSFORM(n_e[*,*,*,it], data,sf,$
                                          NU_SILENT=nu_silent, TOR_FRAC = tor_frac))

  IF (windows_flag) THEN BEGIN ; plot BES locations
      WINDOW, 0
      xr =BES_XLOC[NR_BESchannels/2,0]*Aphys + [-6,6]
      GYRO_RZ_COLOR_CONTOUR, TRANSPOSE(ne_RZ), GYRO_R, GYRO_Z, $
        XRANGE=xr, YRANGE=[-10,2]+BES_Z00, Aphys=Aphys, $
        TITLE = '!4d!Xn!De!N (t = ' + NUMTOSTRING(data.t[it]) + ')'
      
      FOR ir = 0, NR_BESchannels-1 DO FOR iz=0, NZ_BESchannels-1 DO BEGIN
          PLOTS, BES_xloc[ir,iz]*Aphys, BES_zloc[ir,iz]*Aphys, psym=4
          PLOTS, GYRO_R(BES_gyroidx[ir,iz])*Aphys, $
                 GYRO_Z(BES_gyroidx[ir,iz])*Aphys, $
                 psym=2, color=50
      ENDFOR
      CONTOUR, psf.psf, (psf.r+BES_xloc[2,0])*Aphys, $
               (psf.p+BES_zloc[2,0])*Aphys, $
               LEVELS=[0.1,0.5,0.9], $
               C_LINESTYLE=[4,2,0], /OVERPLOT
      CONTOUR, psf_interp[*,*,2,0], GYRO_R*Aphys, GYRO_Z*Aphys, $
               LEVELS=[0.1,0.5,0.9], COLOR=50, $
               C_LINESTYLE=[5,3,1], /OVERPLOT
  ENDIF

  ;set up time axis and result storage arrays
  NT0 = itmax - itmin + 1
  NT = 1 + interp*(NT0-1)
  t = itmin + FINDGEN(NT)/interp
  t_int = itmin+ INDGEN(NT)/interp
  eps = t - t_int

  GYRO_ne_sig = FLTARR(NR_BESchannels, NZ_BESchannels, N_tor_frac, NT)
  syn_ne_sig = FLTARR(NR_BESchannels, NZ_BESchannels, N_tor_frac, NT)
  C_I = COMPLEX(0,1)

  ;begin iterations
  I_THETA = FLTARR(data.n_theta_plot)+1
  FOR it = 0,NT-1 DO BEGIN
      PRINT, it,  '/', NT-1
      ;translate GYRO to rtheta space, apply rho_s renorm
      ne_lf = n_e[*,*,*,t_int[it]]

      ;do linear time interpolation
      IF (eps[it] GT 0) THEN BEGIN
          ne_lf2 = n_e[*,*,*,t_int[it]+1]
          ne_lf = (1. - eps[it])*ne_lf + eps[it]*ne_lf2
      ENDIF

      ;apply doppler shift, and equilibrium density correction
      ;(only does things in global sims)
      FOR i_n = 1, data.n_n-1 DO BEGIN
          ne_lf[*,*,i_n] *= EXP(-C_I*data.n[i_n]*omega0*t[it])/$
                            (REFORM(data.n_eq[data.n_spec-1,*])#I_THETA)
      ENDFOR

      FOR i_tf = 0, n_tor_frac-1 DO BEGIN
          ne_RZ = TRANSPOSE(GYRO_RTHETA_TRANSFORM(ne_lf, data,sf,$
                                                  NU_SILENT=nu_silent, TOR_FRAC=tor_frac))
          
          nu_silent = 1
          FOR ir=0,NR_BESchannels-1 DO FOR iz=0,NZ_BESchannels-1 DO BEGIN
              GYRO_ne_sig[ir,iz,i_tf,it] = ne_RZ[BES_gyroidx[ir,iz]]
              syn_ne_sig[ir,iz,i_tf,it] = TOTAL(ne_RZ*psf_interp[*,*,ir,iz])/$
                psf_interp_norm[ir,iz]
          ENDFOR
      ENDFOR
  ENDFOR

  IF (windows_flag) THEN BEGIN
  	WINDOW, 1
  	PLOT, t, GYRO_ne_sig[0,0,0,*], /XS, title='GYRO ne, syn. BES (red)'
  	OPLOT, t, syn_ne_sig[0,0,0,*], COLOR=100
  ENDIF

  results = {SF:sf, NT:nt, T:t, INTERP:interp, $
                 ITMIN:itmin, ITMAX:itmax, $
                 Aphys:Aphys, omega0:omega0, $
                 N_tor_frac: n_tor_frac, $
                 a_over_cs: a_over_cs, $
                 NR_BESCHANNELS:NR_BESchannels, $
                 NZ_BESCHANNELS:NZ_BESchannels, $
                 GYRO_ne:GYRO_ne_sig, syn_ne:syn_ne_sig, $
                 syn_BES_XLOC:BES_xloc, syn_BES_ZLOC: BES_zloc, $
                 GYRO_BES_XLOC:GYRO_R[BES_gyroidx], $
                 GYRO_BES_ZLOC:GYRO_Z[BES_gyroidx]}

  RETURN, results
END ;make_synarrays
