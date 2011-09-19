FUNCTION read_gyro_profile_data, simdir, profile_data ;, $
  ;THETA_MULT = theta_mult, THETA_PLOT = theta_plot_in
;
; C. Holland, UCSD
;
; v1.0 1-30-2007
;      3-30-2007: added rho_s field
;      5-18-2007: added density, temperature, Er profiles to output
;      struct
;      6-7-2008: added Bunit field
;      6-19-2008: added sim and expt diffusion profiles
;
; Given GYRO sim directory simdir, reads in out.gyro.profile.
; Basically a cut and paste of vugyro routine read_profile_vugyro.pro.
;
; Returns 1 if file read successfull, 0 if not.  Profile data is
; stored in argument profile_data.
;
; v1.1: added theta_mult keyword for reading in nu from
; geometry_arrays.out
;
; v2.: 10-15-2007: added reading in diff_n.out
; v3: 1-3-2008: added theta_plot keyword for use with "reduced" big sims
; v3.1: added diff_to_flow_ei, diff_to_flow_ee, diff_to_flow_ne to
;       output, which are n_r-pt arrays for converting radial
;       chi_i,e,Dne -> q_i,e,gamma_ne and n_bnd for boundary_method=1 sims
; 1-7-2009: updated to use 8-field strucutre of geometry_arrays.out
; 3-3-2009: added version check for nu_geo, compatible with pre and
;           post version 7.0.  Updated for generalized geometry description
; 4-20-2009: added n_moment read
; v4.0: 9-1-2009: compatibility update for GYRO 8.3, generalized ion
; treatment in diffusivity profiles/spectra
; v4.1: 10-26-2009: updated compatibility with geometry_arrays.out, read in flux info
; v5.0: cleaned up to be compatible with GYRO version 9.1.  REmoved
; lots of backwards compatibility, now autmatically gets theta_mult, theta_plot
; 6.15.2010: added V' load from INPUT_profiles.extra
;
  dirpath = GETENV('GYRO_DIR') + '/sim/' + simdir
  filepath =  dirpath + '/out.gyro.profile'  

  ;;get version info
  version_str = ' '
  OPENR, lun, dirpath + '/VERSION_tag', ERR=i_err, /GET_LUN
  IF (i_err EQ 0) THEN BEGIN
      READF, lun, version_str
      FREE_LUN, lun
  ENDIF ELSE BEGIN
      MESSAGE, 'VERSION_tag not found, defaulting to local GYRO', /INFO
      OPENR, lun, GETENV('GYRO_DIR') + '/VERSION', /GET_LUN
      READF, lun, version_str
      FREE_LUN, lun
  ENDELSE

  ;;------------------------------------------------
  ;; Want these set explicitly to integers
  ;;
  n_theta_section = 0
  n_fine	  = 0
  n_r             = 0
  n_kinetic       = 0
  n_pass          = 0
  n_trap          = 0
  n_energy        = 0
  n_theta_plot    = 0
  n_0             = 0
  n_n             = 0
  n_dn            = 0
  nonlinear_flag  = 0
  electron_method = 0
  boundary_method = 0
  ;;------------------------------------------------

  OPENR,1,filepath,ERROR=i_err
  IF (i_err EQ 0) THEN BEGIN
     READF,1,n_r
     READF,1,n_theta_section
     READF,1,n_pass
     READF,1,n_trap
     READF,1,n_energy
     READF,1,n_theta_plot
     READF,1,n_0
     READF,1,n_n
     READF,1,n_dn
     READF,1,n_bnd
     READF,1,nonlinear_flag
     READF,1,electron_method
     READF,1,n_field                   ;1 for phi only, 2 for phi + A_||
     READF,1,n_ion                     ;# ion species
     READF,1,n_kinetic                 ;# evolved species (e- + ion)
     READF,1,n_spec                    ;total # ions and e-
     READF,1,field_r0_flag
     READF,1,field_r0_grid
     READF,1,n_rho
     READF,1,boundary_method

     ;; Need to do some calculations, then allocate 
    ;; next arrays to be read

     n_lambda = n_pass+n_trap

     ;_s indicates variable is surface function, which will be
     ;nonuniform in a non-local simulation, but is held constant
     ;over radial range for flux tube (i.e. "flat profile"?)
     ;ex: in a flux tube T and a/L_T are 2 independant, radially 
     ;uniform constants, but not in a global sim.

     ;n_spec = [ions0, ions1, ... e-]

     r           = fltarr(n_r)          ;midplane minor radius/a
     q         = fltarr(n_r)          ;q(r) even for flux-tube
     r_s         = fltarr(n_r)          ;r, but r0 ~ r[n_r/2] for flat profiles
     q_s         = fltarr(n_r)          ;q, but q0 ~ q[n_r/2] for flat profiles
     dlntdr_s    = fltarr(n_spec,n_r)   ;dlnTdr_s = a/L_T
     dlnndr_s    = fltarr(n_spec,n_r)   ;dlnndr_s = a/L_n
     tem_s       = fltarr(n_spec,n_r)   ;T_s
     den_s       = fltarr(n_spec,n_r)   ;n_s
     phi_doppler_s = fltarr(n_r)          ;equilibrium potential
     aspect_s    = fltarr(n_r)
     delta_s     = fltarr(n_r)
     zeta_s      = fltarr(n_r)
     kappa_s     = fltarr(n_r)
     shift_s     = fltarr(n_r)
     shat_s      = fltarr(n_r)
     s_delta_s   = fltarr(n_r)
     s_zeta_s    = fltarr(n_r)
     s_kappa_s   = fltarr(n_r)
     zmag_s      = fltarr(n_r)
     dzmag_s     = fltarr(n_r)
     beta_unit_s = fltarr(n_r)
     pgamma_s    = fltarr(n_spec,n_r)
     b_unit_s    = fltarr(n_r)
     dr_eodr     = fltarr(n_r)
     z_eff_s     = fltarr(n_r)
     nu_s        = fltarr(n_r)
     gamma_eb_s  = fltarr(n_r)
     w0_s        = fltarr(n_r)

     chi_i_exp       = fltarr(n_r)
     chi_e_exp       = fltarr(n_r)
     diff_to_flow_e1 = fltarr(n_r)
     diff_to_flow_e2 = fltarr(n_r)
     eta_i_tot_exp   = fltarr(n_r)
     aolvi_exp       = fltarr(n_r)
     diff_to_flow_mi = fltarr(n_r) 
     diff_ne_exp     = fltarr(n_r)
     aolne_exp       = fltarr(n_r)
     diff_to_flow_ne = fltarr(n_r)
     eta_i_diff_exp  = fltarr(n_r)
     diff_to_flow_heating = fltarr(n_r)

     lambda  = fltarr(n_lambda)
     energy  = fltarr(n_energy)
     kt_rho  = fltarr(n_n)
     zcharge = fltarr(n_spec)

     ;; basic profile data
     READF,1,r
     READF,1,q
     READF,1,r_s
     READF,1,q_s
     READF,1,dlntdr_s
     READF,1,dlnndr_s
     READF,1,tem_s
     READF,1,den_s
     READF,1,phi_doppler_s   
     READF,1,aspect_s
     READF,1,delta_s
     READF, 1, zeta_s
     READF,1,kappa_s
     READF,1,shift_s
     READF,1,shat_s
     READF,1,s_delta_s
     READF, 1, s_zeta_s
     READF,1,s_kappa_s
     READF, 1, zmag_s
     READF, 1, dzmag_s
     READF,1,beta_unit_s
     READF,1,pgamma_s
     READF,1,b_unit_s
     READF,1,dr_eodr
     READF,1,z_eff_s
     READF,1,nu_s
     READF,1,gamma_eb_s
     READF,1, w0_s    

     readf, 1, box_multiplier ;unused?

     ;; power flows
     READF,1,chi_i_exp
     READF,1,chi_e_exp    
     READF,1,diff_to_flow_e1
     READF,1,diff_to_flow_e2
     READF,1,eta_i_tot_exp
     READF,1,diff_to_flow_mi
     READF,1,aolvi_exp
     READF,1,diff_ne_exp
     READF,1,diff_to_flow_ne
     READF,1,aolne_exp
     READF,1,diff_to_flow_heating

     ;; velocity space, etc
     READF,1,lambda
     READF,1,energy
     READF,1,lambda_tp     
     READF,1,kt_rho
     READF,1,rho_s
     READF,1,zcharge

     READF,1,nfine 
     READF, 1, n_moment
     CLOSE, 1

     theta_mult = nfine/n_theta_plot
;print, nfine, n_theta_plot, theta_mult
;     IF KEYWORD_SET(theta_mult) THEN BEGIN
         geofile = dirpath + '/geometry_arrays.out'
         geo_arr = FLTARR(11, nfine, n_r) ;updated 8->11 10.26.09 per JC update
         exists_nu_geo = READ_GYRO_ARRAY(geo_arr, geofile)
         IF (exists_nu_geo) THEN nu_geo = REFORM(geo_arr[0,*,*]) $
         ELSE nu_geo = FLTARR(nfine,n_r)
;     ENDIF ELSE BEGIN
;         theta_mult = 0.
;         nu_geo = 0.
;         exists_nu_geo = 0
;     ENDELSE

     ; 6-19-2007: read in diffusion profiles
     ; 9-1-2009: add in impurity
     time = READ_GYRO_TIMEVECTOR(dirpath, /SILENT)
     n_time = N_ELEMENTS(time)
     array = FLTARR(n_kinetic,n_field,n_moment,n_r,n_time)
     IF (READ_GYRO_ARRAY(array, dirpath+'/diff_i.out')) THEN BEGIN
         Dn_sim = REFORM(array[*,*,0,*,*],[n_kinetic,n_field,n_r,n_time])
         chi_sim = REFORM(array[*,*,1,*,*],[n_kinetic,n_field,n_r,n_time])
     ENDIF ELSE BEGIN
         PRINT, "Couldn't read " + dirpath + "/diff_i.out"
         Dn_sim = FLTARR(n_kinetic, n_field, n_r, n_time)
         chi_sim = FLTARR(n_kinetic, n_field, n_r, n_time)
     ENDELSE

     ; 10-15-2007 read in diff_n.out- n-resolved, box avg. diffusion
     array = FLTARR(n_kinetic,n_field,n_moment,n_n,n_time)
     IF (READ_GYRO_ARRAY(array, dirpath+'/diff_n.out')) THEN BEGIN
         Dn_kt = REFORM(array[*,*,0,*,*],[n_kinetic,n_field,n_n,n_time])
         chi_kt = REFORM(array[*,*,1,*,*],[n_kinetic,n_field,n_n,n_time])
     ENDIF ELSE BEGIN
         PRINT, "Couldn't read " + dirpath + "/diff_n.out"
         Dn_kt = FLTARR(n_kinetic, n_field, n_n, n_time)
         chi_kt = FLTARR(n_kinetic, n_field, n_n, n_time)
     ENDELSE
      
    ; 10-26-2009 read in fluxes
    array = FLTARR(n_kinetic,n_field,4,n_r,n_time)
    IF (READ_GYRO_ARRAY(array, dirpath+'/gbflux_i.out')) THEN BEGIN
    	 gamma_sim = REFORM(array[*,*,0,*,*],[n_kinetic,n_field,n_r,n_time])
         Q_sim = REFORM(array[*,*,1,*,*],[n_kinetic,n_field,n_r,n_time])
         Pi_sim = REFORM(array[*,*,2,*,*],[n_kinetic,n_field,n_r,n_time])
    ENDIF ELSE BEGIN
         PRINT, "Couldn't read " + dirpath + "/gbflux_i.out"
	 gamma_sim = 0.
	 Q_sim = 0.
	 Pi_sim = 0.
    ENDELSE

    array = FLTARR(n_kinetic,n_field,4,n_n,n_time)
    IF (READ_GYRO_ARRAY(array, dirpath+'/gbflux_n.out')) THEN BEGIN
    	 gamma_n_sim = REFORM(array[*,*,0,*,*],[n_kinetic,n_field,n_n,n_time])
         Q_n_sim = REFORM(array[*,*,1,*,*],[n_kinetic,n_field,n_n,n_time])
         Pi_n_sim = REFORM(array[*,*,2,*,*],[n_kinetic,n_field,n_n,n_time])
    ENDIF ELSE BEGIN
         PRINT, "Couldn't read " + dirpath + "/gbflux_n.out"
	 gamma_n_sim = 0.
	 Q_n_sim = 0.
	 Pi_n_sim = 0.
    ENDELSE
    
    array = FLTARR(13)
    OPENR, 1, dirpath+'/units.out', ERR=i_err
    IF (i_err EQ 0) THEN BEGIN
	FOR ii=0,12 DO BEGIN
		READF, 1, dummy
		array[ii] = dummy
	ENDFOR
	CLOSE, 1
    ENDIF ELSE PRINT, "couldn't find units.out"
    chi_gB = array[8]
    gamma_gB = array[9]
    Q_gB = array[10]
    Pi_gB = array[11]

    ;return useful data in structure, easy to add to as needed
      profile_data = {version:version_str,$  ;simulation version
                     n_r:n_r, $ ;# of radial grid points
                     n_theta_plot:n_theta_plot, $  ;# of _saved_ theta points
                     n_n:n_n, $          ;n = n0 + n_dn*INDGEN(n_n)
                     n_0:n_0, $
                     n_dn:n_dn, $
                     n_field:n_field, $  ;1 for phi only, 2 for phi and A_||
                     n_ion:n_ion, $      ;# ions evolved
                     n_spec:n_spec, $    ;total # ions and e-
                     n_kinetic:n_kinetic, $  ;# of evolved ions and e-
		     n_moment:n_moment, $ ;# of moments saved to diff*.out files
                     ktheta: kt_rho, $     ;k_theta rho_s/i as function of n
                     rho_s: rho_s, $     ;rho_star

                     ;Geometry variables
                     R0:aspect_s*r_s, $  ;R_0(r)/a
                     r:r, $              ;r/a
                     q:q, $            ;q(r) (flux-surface avg?)
                     shear:shat_s, $     ;r/q dq/dr (fs avg or midplane?)
                     kappa:kappa_s, $    ;kappa(r)
                     s_kappa:s_kappa_s, $ ;s_kappa = r dln(kappa)/dr
                     delta:delta_s, $
                     s_delta:s_delta_s, $
                     zeta:zeta_s, $   ;only in 8.0+
                     s_zeta:s_zeta_s, $  ;only in 8.0+
                     zmag:zmag_s, $   ;only in 8.0+
                     dzmag:dzmag_s, $  ;only in 8.0+
                     theta_mult:theta_mult, $
                     nu_geo:nu_geo, $
                     exists_nu_geo:exists_nu_geo, $

                     ;mean profiles
                     n:den_s, $         ;density profiles
                     dlnndr: dlnndr_s, $  ;a/Ln
                     T:tem_s, $         ;temperature profiles
                     dlnTdr: dlnTdr_s, $  ;a/LT
                     Z:zcharge, $       ;Z(species)
                     B_unit:B_unit_s, $  ;effective_B
                     w0:w0_s, $ ;rotation frequency

                     ;transport coeff: note all in gyroBohm units
                     chi_sim: chi_sim, $
                     D_sim: Dn_sim, $
                     chi_kt: chi_kt, $
                     D_kt: Dn_kt, $
                     chi_i_exp: chi_i_exp, $      ;expt. chi_i_profile
                     diff_to_flow_ei: diff_to_flow_e1,$ ;chi_i -> Qi(MW)
                     chi_e_exp: chi_e_exp, $      ;exp. chi_e profile
                     diff_to_flow_ee: diff_to_flow_e2,$ ;chi_e -> Qe(MW)
                     diff_ne_exp: diff_ne_exp, $  ;exp. diff. coeff
                     diff_to_flow_ne: diff_to_flow_ne,$ ;D_ne -> Qi(MW/kev)

                      gamma_sim: gamma_sim, $
                      Q_sim: Q_sim, $
                      Pi_sim: Pi_sim, $
                      gamma_n_sim: gamma_n_sim, $
                      Q_n_sim: Q_n_sim, $
                      Pi_n_sim: Pi_n_sim, $

                      chi_gb:chi_gb, $ ;norm b/w gB (default) and 
                                                  ;m^2/s
                      Gamma_gB: Gamma_gB*624., $ ;10**19/m**2/s
                      Q_gB: Q_gB*100., $ ;W/cm**2
                      Pi_gB: Pi_gB, $ ;NM/m**2

                      n_bnd: n_bnd $ ;# pts in boundary layer
                     }
     ierr = 0
     PRINT, 'Read ' + filepath
   
  ENDIF ELSE BEGIN

     PRINT, 'FATAL ERROR: missing profile_vugyro.out'
     ierr = 1

  ENDELSE

  RETURN, NOT(ierr)
END ;read_gyro_profile_data

