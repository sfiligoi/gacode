FUNCTION read_gyro_profile_data, simdir, profile_data
;
; C. Holland, UCSD
;
; v6.0: 8.25.2011
;         updated for new gacode out.gyro.profile file

  dirpath = GETENV('GYRO_DIR') + '/sim/' + simdir
  filepath =  dirpath + '/out.gyro.profile'  

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
     READF,1,n_field                   ;1 for phi only, 2 for phi + A_||, 3 for phi+A_||+B_||
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
     gamma_e_s = fltarr(n_r)
     gamma_p_s = fltarr(n_r)
     mach_s     = fltarr(n_r)
     b_unit_s    = fltarr(n_r)
     dr_eodr     = fltarr(n_r)
     z_eff_s     = fltarr(n_r)
     nu_s        = fltarr(n_r)
     w0_s        = fltarr(n_r)

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
     READF,1,gamma_e_s
     READF,1,gamma_p_s
     READF,1,mach_s
     READF,1,b_unit_s
     READF,1,dr_eodr
     READF,1,z_eff_s
     READF,1,nu_s
     READF,1, w0_s    

     readf, 1, box_multiplier ;unused?


     ;; velocity space, etc
     READF,1,lambda
     READF,1,energy
     READF,1,lambda_tp     
     READF,1,kt_rho
     READF,1,rho_s
     READF,1,zcharge

     READF,1,nfine 
     CLOSE, 1

     theta_mult = nfine/n_theta_plot
     exists_nu_geo = 0
     nu_geo = 0
     array = FLTARR(11,nfine,n_r)
     OPENR, 1, dirpath+'/out.gyro.geometry_arrays', ERR=i_err
     IF (i_err NE 0) THEN PRINT, "couldn't find out.gyro.geometry_arrays" ELSE BEGIN
         READF, 1, array
         CLOSE, 1
         exists_nu_geo = 1
         nu_geo = REFORM(array[0,*,*])
     ENDELSE 
     
     array = FLTARR(13)
     OPENR, 1, dirpath+'/out.gyro.units', ERR=i_err
     IF (i_err NE 0) THEN PRINT, "couldn't find units.out" ELSE BEGIN
         FOR ii=0,12 DO BEGIN
             READF, 1, dummy
             array[ii] = dummy
         ENDFOR
         CLOSE, 1
     ENDELSE 
     Aphys = array[2]
     csda = array[3]
     Gamma_gB = array[9]
     Q_gB = array[10]
     Pi_gB = array[11]
     S_gB = array[12]

     ;return useful data in structure, easy to add to as needed
     profile_data = {n_r:n_r, $ ;# of radial grid points
                     n_theta_plot:n_theta_plot, $  ;# of _saved_ theta points
                     n_n:n_n, $          ;n = n0 + n_dn*INDGEN(n_n)
                     n_0:n_0, $
                     n_dn:n_dn, $
                     n_field:n_field, $  ;1 for phi only, 2 for phi and A_||
                     n_ion:n_ion, $      ;# ions evolved
                     n_spec:n_spec, $    ;total # ions and e-
                     n_kinetic:n_kinetic, $  ;# of evolved ions and e-
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
                     gamma_e: gamma_e_s, $ ;ExB shear
                     gamma_p: gamma_p_s, $ ;parallel shear
                     mach: mach_s, $;Mach #
                     w0:w0_s, $ ;rotation frequency

                     Aphys: Aphys*100,  $ ;cm
                     csda: csda*1e-3,      $ ;kHz
                     Gamma_gB: Gamma_gB*624., $ ;10**19/m**2/s
                     Q_gB: Q_gB*100., $ ;W/cm**2
                     Pi_gB: Pi_gB, $ ;Nm/m**2
                     S_gB: S_gB, $ ;W/cm**3

                     n_bnd: n_bnd $ ;# pts in boundary layer
                    }
     ierr = 0
     PRINT, 'Read ' + filepath
   
  ENDIF ELSE BEGIN

     PRINT, 'FATAL ERROR: missing profile_vugyro.out'
     ierr = 1

  ENDELSE

  RETURN, ABS(NOT(ierr))
END ;read_gyro_profile_data

