FUNCTION get_gyro_data, simdir, READ_LARGE = read_large, $
  ;THETA_MULT = theta_mult, THETA_PLOT = theta_plot_in, $
  TGYRO=tgyro
;
; C. Holland, UCSD
; v1.0: 1-30-2007
;       3-30-2007: added support for moment_n, moment_e via
;       read_large.  Also added i_field, i_kin flags,i_moment to help specifiy
;       field to analyze
;
; a pared-down version of loadsim.pro from vugyro
; read in profile data (profile_sim.out), time (t.out), phi (u.out, so
; A_|| as well).  Set flags for also reading mom_n.out, mom_e.out.
;
; read_large = 0 (default)- read in u.out, moment_n.out, moment_e.out
;            = 1 u.out only
;            = 2 moment_n.out only
;            = 3 u.out and moment_n.out
;            = 4 moment_n.out and moment_e.out
;            = 5 none
;
; set i_field = 0 to look at phi; = 1 for A (currently unused)
; set i_kin = 0,n_kin-1 to specify which species to look at diff
; species.  e- are i_kin = n_kin-1
;
; v1.1: 4-3-2007: added THETA_MULT keyword to pass to
;       READ_GYRO_PROFILE_DATA for reading in nu from
;       geometry_arrays.out
;
; v2: 10-5-2007: updated read_large to use moment_n and moment_e, 5=out
; v3: 1-3-2008: added theta_plot keyword to change theta_plot for
; reduced big sims
; v4: 4-20-2009: added n_moment functionality
;

  dirpath = GETENV('GYRO_DIR') + '/sim/' + simdir

;  IF (READ_GYRO_PROFILE_DATA(simdir, profile_data, THETA_MULT=theta_mult, $
;                            THETA_PLOT=theta_plot_in)) $
  IF (READ_GYRO_PROFILE_DATA(simdir, profile_data)) THEN BEGIN
      n_r = profile_data.n_r
      n_n = profile_data.n_n
      n_theta_plot = profile_data.n_theta_plot
      
      IF (n_theta_plot EQ 1) THEN theta = 0. ELSE $
        theta = -!PI+2*!PI*FINDGEN(n_theta_plot)/n_theta_plot

      n = profile_data.n_0 + profile_data.n_dn*INDGEN(profile_data.n_n)

      time = READ_GYRO_TIMEVECTOR(dirpath)
      n_time = N_ELEMENTS(time)

      ;read in large field arrays
      DEFAULT, read_large, 0
      exists_phi = 0
      exists_n = 0
      exists_e = 0
      exists_t = 0
      phi = 0
      mom_n = 0
      mom_e = 0
      mom_t = 0

      CASE read_large OF
          0: BEGIN
              exists_phi = 1
              exists_n = 1
              exists_e = 1
          END
          1: BEGIN
              exists_phi = 1
          END
          2: BEGIN
              exists_n = 1
          END
          3: BEGIN
              exists_phi = 1
              exists_n = 1
          END
          4: BEGIN
              exists_n = 1
              exists_e = 1
          END
          5: BREAK
      ENDCASE

      ;start with phi- should fix to include A_||?
      IF (exists_phi) THEN BEGIN
          array = FLTARR(2,n_theta_plot,n_r,profile_data.n_field,n_n,n_time)
          exists_phi = READ_GYRO_ARRAY(array, dirpath+'/u.out')
          IF (exists_phi) THEN BEGIN
              phi = COMPLEXARR(n_theta_plot, n_r, n_n, n_time)
              form = [n_theta_plot, n_r, n_n, n_time]
              phi = COMPLEX(REFORM(TEMPORARY(array[0,*,*,0,*,*]),form),$
                            REFORM(TEMPORARY(array[1,*,*,0,*,*]),form))
          ENDIF
      ENDIF

      form = [n_theta_plot, n_r, profile_data.n_kinetic, n_n, n_time]

      IF (exists_n) THEN BEGIN
          array = FLTARR(2,n_theta_plot,n_r,profile_data.n_kinetic,n_n,n_time)
          exists_n = READ_GYRO_ARRAY(array, dirpath+'/moment_n.out')
          IF (exists_n) THEN BEGIN
              mom_n = COMPLEXARR(n_theta_plot, n_r, profile_data.n_kinetic, $
                                 n_n, n_time)
              mom_n = COMPLEX(REFORM(TEMPORARY(array[0,*,*,*,*,*]),form), $
                              REFORM(TEMPORARY(array[1,*,*,*,*,*]),form))
          ENDIF
      ENDIF

      IF (exists_e) THEN BEGIN
          array = FLTARR(2,n_theta_plot,n_r,profile_data.n_kinetic,n_n,n_time)
          exists_e = READ_GYRO_ARRAY(array, dirpath+'/moment_e.out')
          IF (exists_e) THEN BEGIN
              mom_e = COMPLEXARR(n_theta_plot, n_r, profile_data.n_kinetic, $
                                 n_n, n_time)
              mom_e = COMPLEX(REFORM(TEMPORARY(array[0,*,*,*,*,*]),form), $
                              REFORM(TEMPORARY(array[1,*,*,*,*,*]),form))
          ENDIF    
      ENDIF
      
      ;diffusivity time histories
      array = FLTARR(profile_data.n_kinetic,profile_data.n_field,$
                     profile_data.n_moment,n_time)
      IF (READ_GYRO_ARRAY(array, dirpath+'/diff.out')) THEN BEGIN
          D = FLTARR(profile_data.n_kinetic,profile_data.n_field,n_time)
          chi = FLTARR(profile_data.n_kinetic,profile_data.n_field,n_time)
          D[*,*,*] = array[*,*,0,*]
          chi[*,*,*] = array[*,*,1,*]
      ENDIF ELSE BEGIN
          D = 0
          chi = 0
      ENDELSE

      ;flux time histories
      array = FLTARR(profile_data.n_kinetic,profile_data.n_field,$
                     4,n_time)
      IF (READ_GYRO_ARRAY(array, dirpath+'/gbflux.out')) THEN BEGIN
         Q = FLTARR(profile_data.n_kinetic,profile_data.n_field,n_time)
         Q[*,*,*] = array[*,*,1,*]
      ENDIF ELSE BEGIN
         Q = 0
      ENDELSE
      ;set labels here
      kin_tags = STRARR(profile_data.n_kinetic)
      IF (profile_data.n_kinetic EQ profile_data.n_ion) THEN BEGIN
          ;only have ions
          IF (profile_data.n_kinetic EQ 1) THEN kin_tags[0] = 'ion' $
          ELSE FOR ii = 0, profile_data.n_kinetic-1 DO $
                   kin_tags[ii] = 'ion ' + NUMTOSTRING(ii)
      ENDIF ELSE BEGIN
          ;have ions and electrons
          IF (profile_data.n_kinetic EQ 2) THEN kin_tags[0] = 'ion' $
          ELSE FOR ii = 0, profile_data.n_kinetic-2 DO $
                   kin_tags[ii] = 'ion ' + NUMTOSTRING(ii)
          kin_tags[profile_data.n_kinetic-1] = 'electron'
      ENDELSE

      ;calculate temp fluctuations normalized to Te
      if ((exists_n EQ 1) AND (exists_e EQ 1)) THEN exists_t = 1
      IF (exists_t) THEN BEGIN
          PRINT, 'adding temp flucs'
          mom_t = COMPLEXARR(n_theta_plot, n_r, profile_data.n_kinetic, $
                             n_n, n_time)
          FOR i_r = 0, n_r-1 DO FOR i_spec = 0, profile_data.n_kinetic-1 DO $
            mom_t[*,i_r,i_spec,*,*] = ((2./3)*mom_e[*,i_r,i_spec,*,*] - $
              profile_data.T[i_spec,i_r]*mom_n[*,i_r,i_spec,*,*])/$
              profile_data.n[i_spec, i_r]
      ENDIF
      file = dirpath+'/units.out'

      openr,1,file,err=i_err
      if (i_err eq 0) then begin
        xunits = fltarr(13)
        for i=0,12 do begin
           readf,1,dummy
           xunits[i] = dummy
        endfor
      endif else xunits=0
      close,1
      ;create struct
      gyro_data = {n_r:n_r, n_theta_plot:n_theta_plot, n_n:n_n, n_time:n_time,$
                   r:profile_data.r, theta:theta, n:n, t:time, $
                   profile_data:profile_data, exists_phi:exists_phi, $
                   exists_n: exists_n, exists_e:exists_e, $
                   exists_t: exists_t, $
                   phi:TEMPORARY(phi), mom_n:TEMPORARY(mom_n), $
                   mom_e:TEMPORARY(mom_e), mom_t:TEMPORARY(mom_t), $
                   D:D, chi:chi, Q:Q, kin_tags:kin_tags, $
                   i_field:0, i_kin:0, simdir:simdir, xunits:xunits}     

      RETURN, gyro_data
  ENDIF ELSE BEGIN
      PRINT, "Couldn't read out.gyro.profile, returning 0"
      RETURN, 0
  ENDELSE
END ;get_gyro_data
