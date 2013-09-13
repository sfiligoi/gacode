      SUBROUTINE tglf_input_defaults
!-------------------------------------------------------------------------
! -- set defaults for the tglf model inputs
!----------------------------------------------HSJ-1/2/13-------
    USE nrtype,          ONLY : DP

    USE tglfin,          ONLY :                                                & 
                            adi_elec_tg,adiabatic_elec_tg,                     &
                            find_width_tg,alpha_n_tg,                          &
                            nbasis_max_tg,nb_max_tg,                           &
                            nbasis_min_tg,nb_min_tg,                           &
                            nx_tg,nxgrid_tg,                                   &
                            new_eikonal_tg,ibranch_tg,                         &
                            nmodes_tg, iflux_tg,                               &
                            xky0_tg,ky_tg,                                     &
                            width_max_tg, width_min_tg, nwidth_tg,             &
                            park_tg, ghat_tg, gchat_tg, alpha_e_tg,            &
                            gamma_e_tg, alpha_p_tg, gamma_p_tg, igeo_tg,       &
                            theta_trap_tg,theta_trapped_tg,                    &
                            cnorm_tg, theta0_tg,rmin_tg,                       &
                            rmaj_tg,use_bisection_tg,mass_tg, zs_tg, taus_tg,  &
                            as_tg, rlns_tg, rlts_tg,q_tg, xnuei_tg, wd_zero_tg,&
                            betae_tg, shat_tg, alpha_tg,xwell_tg, kappa_tg,    &
                            s_kappa_tg, delta_tg, s_delta_tg, zeff_tg,         &
                            debye_tg,drmindx_tg,drmajdx_tg,j_surface_tg,       &
                            debye_fac_tg,debye_factor_tg,                      &
                            use_bper_tg,alpha_t_tg,                            &
                            use_bpar_tg,q_prime_tg,                            & 
                            p_prime_tg, filter_tg,                             &
                            Linsker_tg, Linsker_factor_tg,                     &
                            gradB_tg,gradB_factor_tg,                          &
                            nky_tg, b_model_tg, ft_model_tg, xnuei_fac_tg,     &
                            shift_tg,sat_rule_tg,                              &
                            etg_fac_tg,etg_factor_tg,                          &
                            taui_tg,xnu_tg,rlne_tg,                            &
                            rlni_tg,rlnimp_tg,rlte_tg,rlti_tg,aiwt_tg,         &
                            aimpwt_tg,mass_ion_tg,mass_imp_tg,zion_tg,         &
                            zimp_tg,use_mhd_rule_tg,kx0_tg,                    &
                            xnu_fac_tg,xnu_factor_tg,xalpha,                   &
                            shat_mhd_tg,alpha_mhd_tg,                          &
                            x_psi_tg,ialphastab,test_tglf,alpha_dia,           &
                            tglf_implicit_grad,tglf_flux_mult,                 &
                            nspecies_tg,kygrid_model_tg,xnu_model_tg,          &
                            vpar_mach_tg,vexb_mach_tg,tglf_grid_point_parallel,&
                            tglf_wave_no_parallel,alpha_quench_tg,             &
                            tglf_version_no,dzmajdx_tg,ns_tg,use_TM_tg,        &
                            xnue_tg,zmaj_tg,vpar_tg,alpha_kx_e_tg,             &
                            alpha_kx_p_tg,alpha_kx_n_tg,                       &
                            alpha_kx_t_tg,damp_psi_tg,damp_sig_tg,             &
                            vexb_shear_tg, vpar_shear_tg,vns_shear_tg,         &
                            vts_shear_tg, nfourier_tg,fourier_tg,vexb_tg,      &
                            vpar_shear_model_tg,vpar_model_tg,sign_Bt_tg,      &
                            sign_It_tg,s_zeta_tg,zeta_tg,tglf_rmin_sa_in,      &
                            tglf_q_sa_in,tglf_shat_sa_in,tglf_alpha_sa_in,     &
                            tglf_xwell_sa_in,tglf_theta0_sa_in,                &
                            tglf_b_model_sa_in,tglf_ft_model_sa_in 

    USE solcon_gcnmp,ONLY : use_tglf,itran_max
    
    USE common_constants, ONLY : izero,zeroc

      IMPLICIT NONE
      
 
      !-----------------------------------------------------------------------
      ! NOTE on v1.94 aliases: IF defaults on aliases are changed logic at end of
      ! gcnmp_namelist.f90 must be adjusted !!!!!!!!!!!!!!
      !-----------------------------------------------------------------------



      Linsker_tg         = zeroc                ! multiplies Linsker terms (parallel gradient of FLR terms)   
      Linsker_factor_tg  = zeroc                ! v1.93 alias   

      adi_elec_tg        = .FALSE.              ! if .TRUE. use adiabatic electron approximation     
      adiabatic_elec_tg  = .FALSE.              ! v1.93 alias   

 
      alpha_dia          =  1.0_DP              ! multiplier on diamagnetic term of vexb
      alpha_quench_tg    = zeroc                ! ?
      alpha_n_tg         = zeroc                ! ?
      alpha_t_tg         = zeroc                ! ?
      alpha_kx_e_tg      = zeroc                ! ?
      alpha_kx_p_tg      = zeroc                ! ?
      alpha_kx_n_tg      = zeroc                ! ?      
      alpha_kx_t_tg      = zeroc


      aiwt_tg            = 1.0_DP
      aimpwt_tg          = zeroc
      alpha_e_tg         = 1._DP                ! multiplies the ExB velocity shear              
      alpha_mhd_tg       = zeroc                ! Miller geom
      alpha_p_tg         = zeroc                ! multiplies the parallel velocity shear           
      alpha_tg           = zeroc                ! normalized pressure gradient             
      as_tg(:)           = 1.0_DP               ! density N/N0 for each species           
      b_model_tg         = zeroc                ! do 1 or don't 0 use B(theta) factor in the s-alpha geometry k_per**2   
      betae_tg           = zeroc                ! electron beta              
      cnorm_tg           = 1.0_DP

      damp_psi_tg        = zeroc                ! ?
      damp_sig_tg        = zeroc                ! ?

      debye_fac_tg       = 1.0_DP               ! multiplies the debye length      
      debye_factor_tg    = 1.0_DP               ! v1.93 alias (see note above on aliases)  
                                        
      debye_tg           = zeroc                ! ratio of Debye length to the gyro-radius scale rhos0       
      delta_tg           = zeroc                ! Miller geom,triangularity of flux surface            
 
      drmindx_tg         = 1.0_DP               ! Miller geom input
      drmajdx_tg         = zeroc                ! Miller geom input
      dzmajdx_tg         = zeroc                ! Miller geom input

      etg_fac_tg         = 4._DP                ! exponent parameter for ETG saturation rule   
      etg_factor_tg      = 4._DP                ! v1.93 alias (may use 1.25 default?)
     
      filter_tg          = 2.0_DP               ! filter to remove spurious high-frequency magnetic fluctuations       
      find_width_tg      = .TRUE.               ! find the width that maximizes the growthrate         
      ft_model_tg        = 1                    ! trapped fraction model for s-alpha geometry       
      fourier_tg(:,:)    = zeroc
    
      gamma_e_tg         = zeroc                ! ExB velocity shear             
      gamma_p_tg         = zeroc                ! parallel velocity shear             
      gchat_tg           = 1.0_DP               ! multiplies the toroidal drift closure terms        
      ghat_tg            = 1.0_DP               ! multiplies the toroidal drift terms (not closure terms)        
      gradB_tg           = zeroc                ! multiplies parallel gradient of Log(B) terms     
      gradB_factor_tg    = zeroc                ! v1.93 alias (see note above on aliases)  
      ialphastab         = 1                    ! added by hsj
      ibranch_tg         = -1                   ! if -1 order linear mode spectrum by growthrate,        
      iflux_tg           = .TRUE.               ! if .TRUE. compute quasilinear fluxes 
      igeo_tg            =  izero               ! selects s alpha (=0) or Miller(=1) geometry or Fourier rep  (igeo =2) geometry
      j_surface_tg       = izero                ! ?

      kappa_tg           = 1._DP                ! Miller geom,elongation of flux surface   
      kx0_tg             = zeroc                ! ?

      kygrid_model_tg    = 1                    ! select version of ky-grid to use 
      mass_tg(1)         = 2.723E-4_DP          ! mass/m0 for each species me/md
      mass_tg(2)         = 1.0_DP               ! mass/m0 for each species md/md
      mass_tg(3)         = 6.0_DP               ! mass/m0 for each species mc/md
      mass_ion_tg        =  2.0_DP              ! scalar
      mass_imp_tg        =  12.0_DP             ! scalar
         


      nb_max_tg          = 4                    ! maximum number of poloidal hermite basis functions   
      nbasis_max_tg      = 4                    ! v1.93 alias
      nb_min_tg          = 1                    ! minumum number of poloidal hermite basis functions     
      nbasis_min_tg      = 1                    ! v1.93 alias (see note above on aliases) 
    

      new_eikonal_tg     = .TRUE.               ! if .TRUE. compute the eikonal, if .FALSE. re-use eikonal from previous call   

      nfourier_tg        = 4                    ! ?

      nky_tg             = 15                   ! number of modes for high-k sector          
      nmodes_tg          = 2                    ! number of ustable modes to use in computing fluxes (max=4) 
     
      nspecies_tg        = 2                    ! number of species electrons + ions limited to (max=3) 
      ns_tg              = 2                    ! v1.93 alias (see note above on aliases) 

      nwidth_tg          = 21                   ! maximum number of widths to use in search for maximum growthrate    

 
      nx_tg              = 16                   ! number of nodes for the guass-hermite quadrature       
      nxgrid_tg          = 16                   !v1.93  alias (see note above on aliases) 

 
      p_prime_tg         = zeroc                ! pressure derivative wrt poloidal flux           
      park_tg            = 1.0_DP               ! multiplies the parallel gradient operator           
      q_prime_tg         = zeroc                ! safety factor derivative wrt poloidal flux          
      q_tg               = 2.0_DP               ! safety factor     
          
      rlne_tg            = 1.0_DP
      rlni_tg            = 1.0_DP
      rlnimp_tg          = 1.0_DP
      rlns_tg(:)         = 1.0_DP               ! -a0*dLog(N)/drmin for each species                rlte_tg            = 3.0_DP
      rlti_tg            = 3.0_DP
      rlts_tg(:)         = 3.0_DP               ! -a0*dLog(T)/drmin for each species            
      rmaj_tg            = 3.0_DP               ! major radius of local flux surface/a          
      rmin_tg            = 0.5_DP               ! minor radius of local flux surface          
      sat_rule_tg        = 0                    ! 

      s_delta_tg         = zeroc                ! shear in triagularity (Miller input)       
      sign_Bt_tg         = -1.e6_DP             ! MUST be set from statefile or namelist
      sign_It_tg         = -1.e6_DP             ! 1.e6 value is used only to check if set correctly

      s_kappa_tg         = zeroc                ! Miller geom,shear in elongation             
      shat_tg            = 1.0_DP               ! shifted circle geom,magnetic shear                shat_mhd_tg        = 1.0_DP              ! Miller geom,
      shift_tg           = zeroc                ! Miller geom, dRmajor/drminor                      taui_tg            = 1.0_DP              ! THEY ARE PART OF THE NAMELIST DEFAULTS

      s_zeta_tg          = zeroc                ! Miller geom

      taus_tg(:)         = 1.0_DP               ! temperature T/T0 for each species        

      test_tglf          = .FALSE.     ! USE THIS SWITCH TO JUST CALL TGLF WITH AN INPUT FILES, TGLFIN,
                                       ! set up,with above defaults,for a single point call to tglfft_model
                                       ! no transport is done in this case but usefull tglf info is
                                       ! printed out to file test_tglf_output.

      tglf_flux_mult           = 1.0_DP
      tglf_implicit_grad       = .TRUE.
      tglf_grid_point_parallel = .TRUE.  ! parallelize over grid points
      tglf_wave_no_parallel    = .FALSE. ! parallelize over wave vectors (if true)
      tglf_version_no          = 1.94    ! previous version was 1.83,do not change this
                                         ! because  older version routines were removed!

      !the command tglf -h (on Drop) returns this as version number;
      !==========================================================================
      !TGLF r4-232-gf51540 DROP Linux x86_64 ==> call this verion 1.94 in gcnmp
      !==========================================================================


      ! the following are used in put_Fourier_geometry:
      ! they may be aliases for rmin_tg,q_tg,,shat_tg,alpha_tg,xwell_tg,theta0_tg
      ! b_model_tg,ft_model_tg
      tglf_rmin_sa_in          = 0.5_DP
      tglf_q_sa_in             = 2.0_DP
      tglf_shat_sa_in          = 1.0_DP
      tglf_alpha_sa_in         = zeroc
      tglf_xwell_sa_in         = zeroc
      tglf_theta0_sa_in        = zeroc
      tglf_b_model_sa_in       = 1
      tglf_ft_model_sa_in      = 1 

      theta0_tg          = zeroc                ! kx/(shat ky) radial mode number         
      theta_trap_tg      = 0.7                  ! parameter for trapped fraction model           
      theta_trapped_tg   = 0.7                  ! v1.93 alias (see note above on aliases) 

      use_bisection_tg   = .FALSE.              ! use bisection search method to find width that maximizes growthrate       
      use_bpar_tg        = .FALSE.              ! include parallel magnetic field fluctuations           
      use_bper_tg        = .FALSE.              ! include perpendicular magnetic field fluctuations    

      use_mhd_rule_tg    = .TRUE.               ! ?
     
      use_tglf(:)        = izero                ! default tglf to off for all dep. var. 

      use_TM_tg          = .TRUE.               ! not used in gcnmp ?

      vexb_mach_tg       = zeroc                ! ExB velocity mach number               
      vexb_shear_tg      = zeroc                ! ?
      vns_shear_tg(:)    = zeroc                ! ?
      vts_shear_tg(:)    = zeroc                ! ?
      vexb_tg            = zeroc                ! ?
      vpar_shear_tg(:)   = zeroc                ! parallel velocity shear          
      vpar_shear_model_tg = izero               ! ?
      vpar_mach_tg       = zeroc                ! parallel velocity mach number
      vpar_model_tg      = izero                ! ?
      vpar_tg(:)         = zeroc                ! ?

      wd_zero_tg         = 0.1_DP               ! cutoff for curvature drift eigenvalues ABS(wd)>=wd_zero       
      width_max_tg       = 1.65_DP              ! maximum width to use for the search for the maximum growthrate
      width_min_tg       = 0.3_DP               ! minimum width to use for the search for the maximum growthrate

      xalpha             = 1.0_DP               ! s-alpha geometry  stability multiplier
      xnu_tg             =  zeroc               !  ??
      xky0_tg            = 0.3_DP               ! k_theta*rhos for a single k linear stability call to tglf     
      ky_tg              = 0.3_DP               ! v1.93 alias (see note above on aliases)  

      xnu_fac_tg         =  1.0_DP
      xnu_factor_tg      =  1.0_DP              ! v1.93 alias (see note above on aliases) 


      xnu_model_tg       = 2                    ! select version of trapped-passing boundary ei-collision model to use       
      xnue_tg            = zeroc                ! ?

      xnuei_fac_tg       = 1.0_DP               ! multiplies the trapped-passing boundary 
                                                ! terms in the electron-ion collisions  
      xnuei_tg           = zeroc                ! electron-ion collision frequency
      x_psi_tg           = zeroc                ! shift polorization current due to magnetic fluctuations       
      xwell_tg           = zeroc                ! magnetic well parameter             
      zeff_tg            = 1.0_DP               ! effective charge of ions
      zeta_tg            = zeroc                ! Miller geom
      zion_tg            =  1.0_DP              !
      zimp_tg            =  6.0_DP              !
      zs_tg(:)           = 1.0_DP               ! charge/e for each species 
      zmaj_tg            = zeroc                ! used in Miller geometry


   


      RETURN

      END      SUBROUTINE tglf_input_defaults
    

