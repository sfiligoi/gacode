FUNCTION get_gyro_data, simdir, READ_LARGE = read_large, HDF5=hdf5
;
; C. Holland, UCSD
; v5.0: 8.25.2011
;         Brand new version for use with gacode, new out.gyro.xxx file
;         naming convetions
; v5.1: 11.8.2011
;	  Added in support for HDF5 output
;
; SUMMARY: works like loading a simulation into vugyro, but now puts
;         data into an IDL data structure.  Usage is:
;         IDL> data = GET_GYRO_DATA('<simdir>', READ_LARGE=0-5)
;
; INPUTS:
;         simdir: name of GYRO simulation in $GYRO_DIR/sim directory
;         READ_LARGE: option to control which fluctuation fields are loaded
;                    set = 0 to load no fluctuations (default)
;                          = 1 to load all (out.gyro.moment_u,_n,_e,_v)
;                          = 2 to load only out.gyro.moment_u
;                          = 3 to load only out.gyro.moment_n
;                          = 4 to load only out.gyro.moment_u and _n
;                          = 5 to load only out.gyro.moment_n and _e
;
;                    if out.gyro.moment_n and _e are loaded, will also
;                    generate temperature fluctuations using T = (2/3)p-n
;
;           HDF5: read data from out.gyro.initdata.h5,
;           out.gyro.timedata.h5, rather than ASCII files
;

  dirpath  =GETENV('GYRO_DIR') + '/sim/' + simdir

  use_hdf5 = KEYWORD_SET(hdf5) ;flag for using HDF5
  got_profiles = 0
  IF (use_hdf5 EQ 0) THEN BEGIN
      IF (READ_GYRO_PROFILE_DATA(simdir, profile_data) EQ 0) THEN BEGIN
          PRINT, "Couldn't read out.gyro.profile, returning 0"
          RETURN, 0
      ENDIF ELSE got_profiles = 1
  ENDIF ELSE IF (use_hdf5 EQ 1) THEN BEGIN
      IF (READ_GYRO_PROFILE_DATA_HDF5(simdir, profile_data) EQ 0) THEN BEGIN
          PRINT, "Couldn't read out.gyro.initdata.h5, returning 0"
          RETURN, 0
      ENDIF ELSE got_profiles = 1
  ENDIF

  IF (got_profiles EQ 1) THEN BEGIN
      n_r = profile_data.n_r
      n_n = profile_data.n_n
      n_theta_plot = profile_data.n_theta_plot
      
      IF (n_theta_plot EQ 1) THEN theta = 0. ELSE $
        theta = -!PI+2*!PI*FINDGEN(n_theta_plot)/n_theta_plot

      n = profile_data.n_0 + profile_data.n_dn*INDGEN(profile_data.n_n)

      IF (use_hdf5 EQ 1) THEN BEGIN
	tdata = H5_PARSE(dirpath + '/out.gyro.timedata.h5',/READ_DATA)
	time = tdata.t_current._data
      ENDIF ELSE time = READ_GYRO_TIMEVECTOR(dirpath)
      n_time = N_ELEMENTS(time)

      ;read in potentially large fluctuation field arrays
      DEFAULT, read_large, 0
      exists_u = 0
      exists_phi = 0
      exists_Apar = 0
      exists_Bpar = 0
      exists_n = 0
      exists_e = 0
      exists_t = 0
      exists_v = 0
      phi = 0
      Apar = 0
      Bpar = 0
      mom_n = 0
      mom_e = 0
      mom_t = 0
      mom_v = 0

      CASE read_large OF
          0: BREAK ;no data loaded
          1: BEGIN
              exists_u = 1
              exists_n = 1
              exists_e = 1
              exists_v = 1
          END
          2: BEGIN
              exists_u = 1
          END
          3: BEGIN
              exists_n = 1
          END
          4: BEGIN
              exists_u = 1
              exists_n = 1
          END
          5: BEGIN
              exists_n = 1
              exists_e = 1
          END
      ENDCASE

      ;start with out.gyro.moment_u: phi/A_||/B_||
      IF ((use_hdf5 EQ 0) AND (exists_u)) THEN BEGIN
          array = FLTARR(2,n_theta_plot,n_r,profile_data.n_field,n_n,n_time)
          exists_u = READ_GYRO_ARRAY(array, dirpath+'/out.gyro.moment_u')
          IF (exists_u) THEN BEGIN
              form = [n_theta_plot, n_r, n_n, n_time]

              exists_phi = 1
              phi = COMPLEX(REFORM(TEMPORARY(array[0,*,*,0,*,*]),form),$
                            REFORM(TEMPORARY(array[1,*,*,0,*,*]),form))

              IF (profile_data.n_field GT 1) THEN BEGIN
                  exists_Apar = 1
                  Apar = COMPLEX(REFORM(TEMPORARY(array[0,*,*,1,*,*]),form),$
                                REFORM(TEMPORARY(array[1,*,*,1,*,*]),form))
              ENDIF

              IF (profile_data.n_field GT 2) THEN BEGIN
                  exists_Bpar = 1
                  Bpar = COMPLEX(REFORM(TEMPORARY(array[0,*,*,2,*,*]),form),$
                                REFORM(TEMPORARY(array[1,*,*,2,*,*]),form))
              ENDIF
          ENDIF
      ENDIF

      ;load density moment fluctuations
      form = [n_theta_plot, n_r, profile_data.n_kinetic, n_n, n_time]
      IF ((use_hdf5 EQ 0) AND (exists_n)) THEN BEGIN
          array = FLTARR(2,n_theta_plot,n_r,profile_data.n_kinetic,n_n,n_time)
          exists_n = READ_GYRO_ARRAY(array, dirpath+'/out.gyro.moment_n')
          IF (exists_n) THEN BEGIN
              mom_n = COMPLEX(REFORM(TEMPORARY(array[0,*,*,*,*,*]),form), $
                              REFORM(TEMPORARY(array[1,*,*,*,*,*]),form))
          ENDIF
      ENDIF

      ;load energy moment fluctuations, form is same as moment_n
      IF ((use_hdf5 EQ 0) AND (exists_e)) THEN BEGIN
          array = FLTARR(2,n_theta_plot,n_r,profile_data.n_kinetic,n_n,n_time)
          exists_e = READ_GYRO_ARRAY(array, dirpath+'/out.gyro.moment_e')
          IF (exists_e) THEN BEGIN
              mom_e = COMPLEX(REFORM(TEMPORARY(array[0,*,*,*,*,*]),form), $
                              REFORM(TEMPORARY(array[1,*,*,*,*,*]),form))
          ENDIF    
      ENDIF
      
      ;load v_|| moment fluctuations, form is same as moment_n
      IF ((use_hdf5 EQ 0) AND (exists_v)) THEN BEGIN
          array = FLTARR(2,n_theta_plot,n_r,profile_data.n_kinetic,n_n,n_time)
          exists_v = READ_GYRO_ARRAY(array, dirpath+'/out.gyro.moment_v')
          IF (exists_v) THEN BEGIN
              mom_v = COMPLEX(REFORM(TEMPORARY(array[0,*,*,*,*,*]),form), $
                              REFORM(TEMPORARY(array[1,*,*,*,*,*]),form))
          ENDIF    
      ENDIF

      ;load HDF5 fluctuations
      tskip = 200 ;SHOULD GET FROM H5 files in future work
      IF (use_hdf5 EQ 1) THEN BEGIN

          IF (exists_u) THEN BEGIN
              exists_phi = 1
              phi = COMPLEXARR(n_theta_plot,n_r,n_n,n_time)

              IF (profile_data.n_field GT 1) THEN BEGIN
                  exists_Apar = 1
                  Apar = COMPLEXARR(n_theta_plot,n_r,n_n,n_time)
              ENDIF

              IF (profile_data.n_field GT 1) THEN BEGIN
                  exists_Bpar = 1
                  Bpar = COMPLEXARR(n_theta_plot,n_r,n_n,n_time)
              ENDIF
          ENDIF

          IF (exists_n) THEN $
            mom_n = COMPLEXARR(n_theta_plot,n_r,profile_data.n_kinetic,n_n,n_time)

          IF (exists_e) THEN $
            mom_e = COMPLEXARR(n_theta_plot,n_r,profile_data.n_kinetic,n_n,n_time)

          IF (exists_v) THEN $
            mom_v = COMPLEXARR(n_theta_plot,n_r,profile_data.n_kinetic,n_n,n_time)

          FOR it = 0, n_time-1 DO BEGIN
              tlabel = STRCOMPRESS(STRING(it*tskip),/REMOVE_ALL)
              WHILE (STRLEN(tlabel) LT 5) DO tlabel = '0' + tlabel 
              print, tlabel
              flucfile = 'gyro'+tlabel+'.h5'
              flucdata = H5_PARSE(dirpath+'/'+flucfile, /READ)

              IF (exists_u) THEN BEGIN
                  PRINT, 'Loading HDF5 phi, B fluctuations from ' + flucfile
                  phi[*,*,*,it] = COMPLEX($
                                  TRANSPOSE(flucdata.phi_modes.phi_real._data[*,*,0:n_theta_plot-1],[2,1,0]),$
                                 TRANSPOSE(flucdata.phi_modes.phi_imag._data[*,*,0:n_theta_plot-1],[2,1,0]))

                  IF (exists_Apar) THEN BEGIN
                      Apar[*,*,*,it] = COMPLEX($
                                      TRANSPOSE(flucdata.Apar_modes.Apar_real._data[*,*,0:n_theta_plot-1],[2,1,0]),$
                                      TRANSPOSE(flucdata.Apar_modes.Apar_imag._data[*,*,0:n_theta_plot-1],[2,1,0]))
                  ENDIF

                  IF (exists_Bpar) THEN BEGIN
                      Bpar[*,*,*,it] = COMPLEX($
                                      TRANSPOSE(flucdata.Bpar_modes.Bpar_real._data[*,*,0:n_theta_plot-1],[2,1,0]),$
                                      TRANSPOSE(flucdata.Bpar_modes.Bpar_imag._data[*,*,0:n_theta_plot-1],[2,1,0]))
                  ENDIF
              ENDIF
              
              IF (exists_n) THEN BEGIN
                  PRINT, 'Loading HDF5 density fluctuations from ' + flucfile
                  
                  tmp = COMPLEX($
                        TRANSPOSE(flucdata.density_ion0_modes.density_ion0_real._data,[2,1,0]),$
                        TRANSPOSE(flucdata.density_ion0_modes.density_ion0_imag._data,[2,1,0]))
                  mom_n[*,*,0,*,it] = $
                    REFORM(tmp[0:n_theta_plot-1,*,*], [n_theta_plot,n_r,1,n_n,1])

                  IF (profile_data.n_ion GE 2) THEN BEGIN
                      tmp = COMPLEX($
                            TRANSPOSE(flucdata.density_ion1_modes.density_ion1_real._data,$
                                      [2,1,0]),$
                            TRANSPOSE(flucdata.density_ion1_modes.density_ion1_imag._data,$
                                      [2,1,0]))
                      mom_n[*,*,1,*,it] = $
                        REFORM(tmp[0:n_theta_plot-1,*,*], [n_theta_plot,n_r,1,n_n,1])
                  ENDIF

                  IF (profile_data.n_ion GE 3) THEN BEGIN
                      tmp = COMPLEX($
                            TRANSPOSE(flucdata.density_ion2_modes.density_ion2_real._data, $
                                      [2,1,0]),$
                            TRANSPOSE(flucdata.density_ion2_modes.density_ion2_imag._data, $
                                      [2,1,0]))
                      mom_n[*,*,2,*,it] = $
                        REFORM(tmp[0:n_theta_plot-1,*,*], [n_theta_plot,n_r,1,n_n,1])
                  ENDIF

                  IF (profile_data.n_ion GE 3) THEN BEGIN
                      tmp = COMPLEX($
                            TRANSPOSE(flucdata.density_ion3_modes.density_ion3_real._data, $
                                      [2,1,0]),$
                            TRANSPOSE(flucdata.density_ion3_modes.density_ion3_imag._data,$
                                      [2,1,0]))
                      mom_n[*,*,1,*,it] = $
                        REFORM(tmp[0:n_theta_plot-1,*,*], [n_theta_plot,n_r,1,n_n,1])
                  ENDIF

                  IF (profile_data.n_kinetic GT profile_data.n_ion) THEN BEGIN
                      tmp = COMPLEX(TRANSPOSE(flucdata.density_electron_modes.$
                      density_electron_real._data,[2,1,0]),$
                        TRANSPOSE(flucdata.density_electron_modes.$
                      density_electron_imag._data,[2,1,0]))
                      mom_n[*,*,profile_data.n_kinetic-1,*,it] = $
                        REFORM(tmp[0:n_theta_plot-1,*,*], [n_theta_plot,n_r,1,n_n,1])
                  ENDIF
              ENDIF ;mom_n

              IF (exists_e) THEN BEGIN
                  PRINT, 'Loading HDF5 energy fluctuations from ' + flucfile
                  
                  tmp = COMPLEX($
                        TRANSPOSE(flucdata.energy_ion0_modes.energy_ion0_real._data,[2,1,0]),$
                        TRANSPOSE(flucdata.energy_ion0_modes.energy_ion0_imag._data,[2,1,0]))
                  mom_e[*,*,0,*,it] = $
                    REFORM(tmp[0:n_theta_plot-1,*,*], [n_theta_plot,n_r,1,n_n,1])

                  IF (profile_data.n_ion GE 2) THEN BEGIN
                      tmp = COMPLEX($
                            TRANSPOSE(flucdata.energy_ion1_modes.energy_ion1_real._data,$
                                      [2,1,0]),$
                            TRANSPOSE(flucdata.energy_ion1_modes.energy_ion1_imag._data,$
                                      [2,1,0]))
                      mom_e[*,*,1,*,it] = $
                        REFORM(tmp[0:n_theta_plot-1,*,*], [n_theta_plot,n_r,1,n_n,1])
                  ENDIF

                  IF (profile_data.n_ion GE 3) THEN BEGIN
                      tmp = COMPLEX($
                            TRANSPOSE(flucdata.energy_ion2_modes.energy_ion2_real._data, $
                                      [2,1,0]),$
                            TRANSPOSE(flucdata.energy_ion2_modes.energy_ion2_imag._data, $
                                      [2,1,0]))
                      mom_e[*,*,2,*,it] = $
                        REFORM(tmp[0:n_theta_plot-1,*,*], [n_theta_plot,n_r,1,n_n,1])
                  ENDIF

                  IF (profile_data.n_ion GE 3) THEN BEGIN
                      tmp = COMPLEX($
                            TRANSPOSE(flucdata.energy_ion3_modes.energy_ion3_real._data, $
                                      [2,1,0]),$
                            TRANSPOSE(flucdata.energy_ion3_modes.energy_ion3_imag._data,$
                                      [2,1,0]))
                      mom_e[*,*,1,*,it] = $
                        REFORM(tmp[0:n_theta_plot-1,*,*], [n_theta_plot,n_r,1,n_n,1])
                  ENDIF

                  IF (profile_data.n_kinetic GT profile_data.n_ion) THEN BEGIN
                      tmp = COMPLEX(TRANSPOSE(flucdata.energy_electron_modes.$
                      energy_electron_real._data,[2,1,0]),$
                        TRANSPOSE(flucdata.energy_electron_modes.$
                      energy_electron_imag._data,[2,1,0]))
                      mom_e[*,*,profile_data.n_kinetic-1,*,it] = $
                        REFORM(tmp[0:n_theta_plot-1,*,*], [n_theta_plot,n_r,1,n_n,1])
                  ENDIF
              ENDIF ;mom_e

              IF (exists_v) THEN BEGIN
                  PRINT, 'Loading HDF5 vpar fluctuations from ' + flucfile
                  
                  tmp = COMPLEX($
                        TRANSPOSE(flucdata.v_par_ion0_modes.v_par_ion0_real._data,[2,1,0]),$
                        TRANSPOSE(flucdata.v_par_ion0_modes.v_par_ion0_imag._data,[2,1,0]))
                  mom_v[*,*,0,*,it] = $
                    REFORM(tmp[0:n_theta_plot-1,*,*], [n_theta_plot,n_r,1,n_n,1])

                  IF (profile_data.n_ion GE 2) THEN BEGIN
                      tmp = COMPLEX($
                            TRANSPOSE(flucdata.v_par_ion1_modes.v_par_ion1_real._data,$
                                      [2,1,0]),$
                            TRANSPOSE(flucdata.v_par_ion1_modes.v_par_ion1_imag._data,$
                                      [2,1,0]))
                      mom_v[*,*,1,*,it] = $
                        REFORM(tmp[0:n_theta_plot-1,*,*], [n_theta_plot,n_r,1,n_n,1])
                  ENDIF

                  IF (profile_data.n_ion GE 3) THEN BEGIN
                      tmp = COMPLEX($
                            TRANSPOSE(flucdata.v_par_ion2_modes.v_par_ion2_real._data, $
                                      [2,1,0]),$
                            TRANSPOSE(flucdata.v_par_ion2_modes.v_par_ion2_imag._data, $
                                      [2,1,0]))
                      mom_v[*,*,2,*,it] = $
                        REFORM(tmp[0:n_theta_plot-1,*,*], [n_theta_plot,n_r,1,n_n,1])
                  ENDIF

                  IF (profile_data.n_ion GE 3) THEN BEGIN
                      tmp = COMPLEX($
                            TRANSPOSE(flucdata.v_par_ion3_modes.v_par_ion3_real._data, $
                                      [2,1,0]),$
                            TRANSPOSE(flucdata.v_par_ion3_modes.v_par_ion3_imag._data,$
                                      [2,1,0]))
                      mom_v[*,*,1,*,it] = $
                        REFORM(tmp[0:n_theta_plot-1,*,*], [n_theta_plot,n_r,1,n_n,1])
                  ENDIF

                  IF (profile_data.n_kinetic GT profile_data.n_ion) THEN BEGIN
                      tmp = COMPLEX(TRANSPOSE(flucdata.v_par_electron_modes.$
                      v_par_electron_real._data,[2,1,0]),$
                        TRANSPOSE(flucdata.v_par_electron_modes.$
                      v_par_electron_imag._data,[2,1,0]))
                      mom_v[*,*,profile_data.n_kinetic-1,*,it] = $
                        REFORM(tmp[0:n_theta_plot-1,*,*], [n_theta_plot,n_r,1,n_n,1])
                  ENDIF
              ENDIF ;mom_v

          ENDFOR
      ENDIF ;HDF5 flucs

      ;calculate temp fluctuations normalized to Te
      IF ((exists_n EQ 1) AND (exists_e EQ 1)) THEN exists_t = 1
      IF (exists_t) THEN BEGIN
          PRINT, 'Generating temperature fluctuations'
          mom_t = COMPLEXARR(n_theta_plot, n_r, profile_data.n_kinetic, $
                             n_n, n_time)
          FOR i_r = 0, n_r-1 DO FOR i_spec = 0, profile_data.n_kinetic-1 DO $
            mom_t[*,i_r,i_spec,*,*] = ((2./3)*mom_e[*,i_r,i_spec,*,*] - $
              profile_data.T[i_spec,i_r]*mom_n[*,i_r,i_spec,*,*])/$
              profile_data.n[i_spec, i_r]
      ENDIF
            
      ;set species labels here
      kin_tags = STRARR(profile_data.n_kinetic)
      IF (profile_data.n_kinetic EQ profile_data.n_ion) THEN BEGIN ;only have ions
          IF (profile_data.n_kinetic EQ 1) THEN kin_tags[0] = 'ion' $
          ELSE FOR ii = 0, profile_data.n_kinetic-1 DO $
                   kin_tags[ii] = 'ion ' + NUMTOSTRING(ii)
      ENDIF ELSE BEGIN ;have ions and electrons
          IF (profile_data.n_kinetic EQ 2) THEN kin_tags[0] = 'ion' $
          ELSE FOR ii = 0, profile_data.n_kinetic-2 DO $
                   kin_tags[ii] = 'ion ' + NUMTOSTRING(ii)
          kin_tags[profile_data.n_kinetic-1] = 'electron'
      ENDELSE

      ;load fluxes- initalize first
      Gamma_t = 0
      Q_t = 0
      Pi_t = 0
      Sexch_t = 0
      Gamma_rt = 0
      Q_rt = 0
      Pi_rt = 0
      Sexch_rt = 0
      Gamma_nt = 0
      Q_nt = 0
      Pi_nt = 0
      Sexch_nt = 0
 
     IF ((use_hdf5 EQ 1) AND (profile_data.nonlinear_flag EQ 1)) THEN BEGIN
         ;load total fluxes vs. time
         dsize = SIZE(tdata.gbflux._data)
         refarr = [dsize[1],dsize[3],dsize[4]]
         Gamma_t = TRANSPOSE(REFORM(tdata.gbflux._data[*,0,*,*],refarr),[2,1,0])
         Q_t = TRANSPOSE(REFORM(tdata.gbflux._data[*,1,*,*],refarr),[2,1,0])
         Pi_t = TRANSPOSE(REFORM(tdata.gbflux._data[*,2,*,*],refarr),[2,1,0])
         Secxh_t = TRANSPOSE(REFORM(tdata.gbflux._data[*,3,*,*],refarr),[2,1,0])

         ;load total fluxes vs.radius & time
         dsize = SIZE(tdata.gbflux_i._data)
         refarr = [dsize[1],dsize[2],dsize[4],dsize[5]]
         Gamma_rt = TRANSPOSE(REFORM(tdata.gbflux_i._data[*,*,0,*,*],refarr),[3,2,1,0])
         Q_rt = TRANSPOSE(REFORM(tdata.gbflux_i._data[*,*,1,*,*],refarr),[3,2,1,0])
         Pi_rt = TRANSPOSE(REFORM(tdata.gbflux_i._data[*,*,2,*,*],refarr),[3,2,1,0])
         Sexch_rt = TRANSPOSE(REFORM(tdata.gbflux_i._data[*,*,3,*,*],refarr),[3,2,1,0])
        
         ;load total fluxes vs. modenumber & time
         Gamma_nt = FLTARR(profile_data.n_kinetic,profile_data.n_field,n_n,n_time)
         Q_nt = FLTARR(profile_data.n_kinetic,profile_data.n_field,n_n,n_time)
         Pi_nt = FLTARR(profile_data.n_kinetic,profile_data.n_field,n_n,n_time)
         Sexch_nt = FLTARR(profile_data.n_kinetic,profile_data.n_field,n_n,n_time)         

         ;load electrostatic component of fluxes vs. modenumber & time
         FOR i_kin = 0, profile_data.n_kinetic-1 DO BEGIN
             i_kin_str = NUMTOSTRING(i_kin)
             com = 'Gamma_nt['+i_kin_str+',0,*,*] = TRANSPOSE(tdata.gbflux_n_ion' + i_kin_str + '_phi_density._data)'
             r = EXECUTE(com)
             com = 'Q_nt['+i_kin_str+',0,*,*] = TRANSPOSE(tdata.gbflux_n_ion' + i_kin_str + '_phi_energy._data)'
             r = EXECUTE(com)
             com = 'Pi_nt['+i_kin_str+',0,*,*] = TRANSPOSE(tdata.gbflux_n_ion' + i_kin_str + '_phi_momentum._data)'
             r = EXECUTE(com)
             com = 'Sexch_nt['+i_kin_str+',0,*,*] = TRANSPOSE(tdata.gbflux_n_ion' + i_kin_str + '_phi_energyexchange._data)'
             r = EXECUTE(com)
         ENDFOR

         ;load Apar component of fluxes vs. modenumber & time
         IF (profile_data.n_field GE 2) THEN FOR i_kin = 0, profile_data.n_kinetic-1 DO BEGIN
             i_kin_str = NUMTOSTRING(i_kin)
             com = 'Gamma_nt['+i_kin_str+',1,*,*] = TRANSPOSE(tdata.gbflux_n_ion' + i_kin_str + '_apar_density._data)'
             r = EXECUTE(com)
             com = 'Q_nt['+i_kin_str+',1,*,*] = TRANSPOSE(tdata.gbflux_n_ion' + i_kin_str + '_apar_energy._data)'
             r = EXECUTE(com)
             com = 'Pi_nt['+i_kin_str+',1,*,*] = TRANSPOSE(tdata.gbflux_n_ion' + i_kin_str + '_apar_momentum._data)'
             r = EXECUTE(com)
             com = 'Sexch_nt['+i_kin_str+',1,*,*] = TRANSPOSE(tdata.gbflux_n_ion' + i_kin_str + '_apar_energyexchange._data)'
             r = EXECUTE(com)
         ENDFOR

        ;load Bpar component of fluxes vs. modenumber & time
        IF (profile_data.n_field EQ 3) THEN FOR i_kin = 0, profile_data.n_kinetic-1 DO BEGIN
             i_kin_str = NUMTOSTRING(i_kin)
             com = 'Gamma_nt['+i_kin_str+',2,*,*] = TRANSPOSE(tdata.gbflux_n_ion' + i_kin_str + '_bpar_density._data)'
             r = EXECUTE(com)
             com = 'Q_nt['+i_kin_str+',2,*,*] = TRANSPOSE(tdata.gbflux_n_ion' + i_kin_str + '_bpar_energy._data)'
             r = EXECUTE(com)
             com = 'Pi_nt['+i_kin_str+',2,*,*] = TRANSPOSE(tdata.gbflux_n_ion' + i_kin_str + '_bpar_momentum._data)'
             r = EXECUTE(com)
             com = 'Sexch_nt['+i_kin_str+',2,*,*] = TRANSPOSE(tdata.gbflux_n_ion' + i_kin_str + '_bpar_energyexchange._data)'
             r = EXECUTE(com)
         ENDFOR

     ENDIF ELSE BEGIN           ;from ascii files
          ;load total fluxes vs. time
          array = FLTARR(profile_data.n_kinetic,profile_data.n_field,4,n_time)
          exists_flux = READ_GYRO_ARRAY(array, dirpath+'/out.gyro.gbflux')
          IF (exists_flux) THEN BEGIN
              form = [profile_data.n_kinetic,profile_data.n_field,n_time]
              Gamma_t = REFORM(array[*,*,0,*],form)
              Q_t = REFORM(array[*,*,1,*],form)
              Pi_t = REFORM(array[*,*,2,*],form)
              Sexch_t= REFORM(array[*,*,3,*],form)
          ENDIF

          ;load total fluxes vs. radius & time
          array = FLTARR(profile_data.n_kinetic,profile_data.n_field,4,n_r,n_time)
          exists_flux = READ_GYRO_ARRAY(array, dirpath+'/out.gyro.gbflux_i')
          IF (exists_flux) THEN BEGIN
              form = [profile_data.n_kinetic,profile_data.n_field,n_r,n_time]
              Gamma_rt = REFORM(array[*,*,0,*,*],form)
              Q_rt = REFORM(array[*,*,1,*,*],form)
              Pi_rt = REFORM(array[*,*,2,*,*],form)
              Sexch_rt= REFORM(array[*,*,3,*,*],form)
          ENDIF

          ;load total fluxes vs. modenumbers & time
          array = FLTARR(profile_data.n_kinetic,profile_data.n_field,4,n_n,n_time)
          exists_flux = READ_GYRO_ARRAY(array, dirpath+'/out.gyro.gbflux_n')
          IF (exists_flux) THEN BEGIN
              form = [profile_data.n_kinetic,profile_data.n_field,n_n,n_time]
              Gamma_nt = REFORM(array[*,*,0,*,*],form)
              Q_nt = REFORM(array[*,*,1,*,*],form)
              Pi_nt = REFORM(array[*,*,2,*,*],form)
              Sexch_nt= REFORM(array[*,*,3,*,*],form)
          ENDIF
      ENDELSE
      
      ;create struct
      gyro_data = {n_r:n_r, $
                   n_bnd: profile_data.n_bnd, $ ;# pts in boundary layer
                   n_theta_plot:profile_data.n_theta_plot, $
                   n_n: profile_data.n_n, $
                   n_0: profile_data.n_0, $
                   n_dn: profile_data.n_dn, $
                   n_field: profile_data.n_field, $ ;1 for phi only, 2 for phi and A_||
                   n_ion: profile_data.n_ion, $ ;# ions evolved
                   n_spec: profile_data.n_spec, $ ;total # ions and e-
                   n_kinetic: profile_data.n_kinetic, $ ;# of evolved ions and e-
                   n_time:n_time,$
                   r: profile_data.r, $
                   theta:theta, $
                   n: n, $
                   ktheta: profile_data.ktheta, $ ;k_theta rho_s/i as function of n
                   t: time, $

                   ;Geometry variables
                   R0: profile_data.R0, $ ;R_0(r)/a
                   q: profile_data.q, $       ;q(r) (flux-surface avg)
                   shear: profile_data.shear, $ ;r/q dq/dr (fs avg)
                   kappa: profile_data.kappa, $ ;kappa(r)
                   s_kappa: profile_data.s_kappa, $ ;s_kappa = r dln(kappa)/dr
                   delta: profile_data.delta, $
                   s_delta: profile_data.s_delta, $
                   zeta: profile_data.zeta, $ 
                   s_zeta: profile_data.s_zeta, $ 
                   zmag: profile_data.zmag, $ 
                   dzmag: profile_data.dzmag, $
                   theta_mult: profile_data.theta_mult, $
                   nu_geo: profile_data.nu_geo, $
                   exists_nu_geo: profile_data.exists_nu_geo, $

                   ;mean profiles
                   n_eq: profile_data.n, $   ;density profiles
                   dlnndr:  profile_data.dlnndr, $ ;a/Ln
                   T_eq: profile_data.T, $   ;temperature profiles
                   dlnTdr:  profile_data.dlnTdr, $ ;a/LT
                   Z: profile_data.Z, $ ;Z(species)
                   B_unit: profile_data.B_unit, $ ;effective_B
                   gamma_e:  profile_data.gamma_e, $ ;ExB shear
                   gamma_p:  profile_data.gamma_p, $ ;parallel shear
                   mach:  profile_data.mach, $ ;Mach #
                   w0: profile_data.w0, $   ;rotation frequency

                   ;normalizations
                   Aphys: profile_data.Aphys, $ ;a in cm
                   csda: profile_data.csda,     $ ;c_s/a in kHz
                   rho_s: profile_data.rho_s, $ ;rho_star
                   Gamma_gB:  profile_data.Gamma_gB, $ ;10**19/m**2/s
                   Q_gB:  profile_data.Q_gB, $ ;W/cm**2
                   Pi_gB:  profile_data.Pi_gB, $ ;Nm/m**2
                   S_gB:  profile_data.S_gB, $ ;W/cm**3

                   ;fluxes
                   Gamma_t: Gamma_t, $ ;total fluxes vs time in gB units
                   Q_t: Q_t, $
                   Pi_t: Pi_t, $
                   Sexch_t: Sexch_t, $
                   Gamma_rt: Gamma_rt, $ ;fluxes vs radius & time in gB units
                   Q_rt: Q_rt, $
                   Pi_rt: Pi_rt, $
                   Sexch_rt: Sexch_rt, $
                   Gamma_nt: Gamma_nt, $ ;fluxes vs modenumber & time in gB units
                   Q_nt: Q_nt, $
                   Pi_nt: Pi_nt, $
                   Sexch_nt: Sexch_nt, $

                   ;fluctuation fields
                   exists_phi:exists_phi, $
                   exists_Apar:exists_Apar, $
                   exists_Bpar:exists_Bpar, $
                   exists_n: exists_n, $
                   exists_e:exists_e, $
                   exists_t: exists_t, $
                   exists_v: exists_v, $
                   phi:TEMPORARY(phi), $
                   Apar:TEMPORARY(Apar), $
                   Bpar:TEMPORARY(Bpar), $
                   mom_n:TEMPORARY(mom_n), $
                   mom_e:TEMPORARY(mom_e), $
                   mom_t:TEMPORARY(mom_t), $
                   mom_v:TEMPORARY(mom_v), $
                   kin_tags:kin_tags, i_field:0, i_kin:0, simdir:simdir}     

      RETURN, gyro_data
  ENDIF
END ;get_gyro_data

