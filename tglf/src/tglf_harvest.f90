  ! Harvest LOCAL INTERFACE variables
  SUBROUTINE tglf_harvest_local()

    USE tglf_interface

    INTEGER :: ierr, i
    CHARACTER(LEN=65507) :: harvest_sendline
    CHARACTER(LEN=2) :: NUM
    CHARACTER NUL
    PARAMETER(NUL = CHAR(0))

    IF (.NOT.tglf_use_transport_model_in) THEN
        RETURN
    ENDIF

    IF (.NOT.tglf_new_eikonal_in) THEN
        RETURN
    ENDIF

    IF (tglf_vexb_in.NE.0) THEN
        RETURN
    ENDIF

    IF (.NOT.tglf_iflux_in) THEN
        RETURN
    ENDIF

    IF (tglf_geometry_flag_in .NE. 1 ) THEN
       WRITE(1,*) 'HARVEST ONLY SUPPORTS MILLER GEOMETRY'
       CLOSE(1)
       RETURN
    ENDIF

    ierr=set_harvest_verbose(1)
    ierr=set_harvest_table('TGLF_harvest?'//NUL)
    ierr=set_harvest_host('127.0.0.1'//NUL)
    ierr=set_harvest_port(32000)
    ierr=init_harvest(harvest_sendline,LEN(harvest_sendline))

!   '#---------------------------------------------------'
!   '# Plasma parameters:'
!   '#---------------------------------------------------'
    ierr=set_harvest_payload_bol(harvest_sendline,'SIGN_BT'//NUL,INT((tglf_sign_bt_in+1)/2.))
    ierr=set_harvest_payload_bol(harvest_sendline,'SIGN_IT'//NUL,INT((tglf_sign_it_in+1)/2.))
    ierr=set_harvest_payload_dbl(harvest_sendline,'KY'//NUL,tglf_ky_in)
    ierr=set_harvest_payload_dbl(harvest_sendline,'VEXB_SHEAR'//NUL,tglf_vexb_shear_in)
    ierr=set_harvest_payload_dbl(harvest_sendline,'BETAE'//NUL,tglf_betae_in)
    ierr=set_harvest_payload_dbl(harvest_sendline,'XNUE'//NUL,tglf_xnue_in)
    ierr=set_harvest_payload_dbl(harvest_sendline,'ZEFF'//NUL,tglf_zeff_in)
    ierr=set_harvest_payload_dbl(harvest_sendline,'DEBYE'//NUL,tglf_debye_in)

!   '#---------------------------------------------------'
!   '# Species vectors:'
!   '#---------------------------------------------------'
    DO i = 1,tglf_ns_in
      IF (i < 10) THEN
         write (NUM, "(I01,A1)") i,NUL
      else
         write (NUM, "(I02,A1)") i,NUL
      ENDIF
      ierr=set_harvest_payload_dbl(harvest_sendline,'ZS_'//NUM,tglf_zs_in(i))
      ierr=set_harvest_payload_dbl(harvest_sendline,'MASS_'//NUM,tglf_mass_in(i))
      ierr=set_harvest_payload_dbl(harvest_sendline,'RLNS_'//NUM,tglf_rlns_in(i))
      ierr=set_harvest_payload_dbl(harvest_sendline,'RLTS_'//NUM,tglf_rlts_in(i))
      ierr=set_harvest_payload_dbl(harvest_sendline,'TAUS_'//NUM,tglf_taus_in(i))
      ierr=set_harvest_payload_dbl(harvest_sendline,'AS_'//NUM,tglf_as_in(i))
      ierr=set_harvest_payload_dbl(harvest_sendline,'VPAR_'//NUM,tglf_vpar_in(i))
      ierr=set_harvest_payload_dbl(harvest_sendline,'VPAR_SHEAR_'//NUM,tglf_vpar_shear_in(i))
      ierr=set_harvest_payload_dbl(harvest_sendline,'VNS_SHEAR_'//NUM,tglf_vns_shear_in(i))
      ierr=set_harvest_payload_dbl(harvest_sendline,'VTS_SHEAR_'//NUM,tglf_vts_shear_in(i))
    ENDDO

!   '#---------------------------------------------------'
!   '# Miller geometry parameters:'
!   '#---------------------------------------------------'
    ierr=set_harvest_payload_dbl(harvest_sendline,'RMIN_LOC'//NUL,tglf_rmin_loc_in)
    ierr=set_harvest_payload_dbl(harvest_sendline,'RMAJ_LOC'//NUL,tglf_rmaj_loc_in)
    ierr=set_harvest_payload_dbl(harvest_sendline,'ZMAJ_LOC'//NUL,tglf_zmaj_loc_in)
    ierr=set_harvest_payload_dbl(harvest_sendline,'DRMINDX_LOC'//NUL,tglf_drmindx_loc_in)
    ierr=set_harvest_payload_dbl(harvest_sendline,'DRMAJDX_LOC'//NUL,tglf_drmajdx_loc_in)
    ierr=set_harvest_payload_dbl(harvest_sendline,'DZMAJDX_LOC'//NUL,tglf_dzmajdx_loc_in)
    ierr=set_harvest_payload_dbl(harvest_sendline,'Q_LOC'//NUL,tglf_q_loc_in)
    ierr=set_harvest_payload_dbl(harvest_sendline,'KAPPA_LOC'//NUL,tglf_kappa_loc_in)
    ierr=set_harvest_payload_dbl(harvest_sendline,'S_KAPPA_LOC'//NUL,tglf_s_kappa_loc_in)
    ierr=set_harvest_payload_dbl(harvest_sendline,'DELTA_LOC'//NUL,tglf_delta_loc_in)
    ierr=set_harvest_payload_dbl(harvest_sendline,'S_DELTA_LOC'//NUL,tglf_s_delta_loc_in)
    ierr=set_harvest_payload_dbl(harvest_sendline,'ZETA_LOC'//NUL,tglf_zeta_loc_in)
    ierr=set_harvest_payload_dbl(harvest_sendline,'S_ZETA_LOC'//NUL,tglf_s_zeta_loc_in)
    ierr=set_harvest_payload_dbl(harvest_sendline,'P_PRIME_LOC'//NUL,tglf_p_prime_loc_in)
    ierr=set_harvest_payload_dbl(harvest_sendline,'Q_PRIME_LOC'//NUL,tglf_q_prime_loc_in)
    ierr=set_harvest_payload_dbl(harvest_sendline,'KX0_LOC'//NUL,tglf_kx0_loc_in)

!   '#---------------------------------------------------'
!   '# Gaussian width parameters:'
!   '#---------------------------------------------------'
    ierr=set_harvest_payload_dbl(harvest_sendline,'+WIDTH'//NUL,tglf_width_in)
    ierr=set_harvest_payload_dbl(harvest_sendline,'+WIDTH_MIN'//NUL,tglf_width_min_in)
    ierr=set_harvest_payload_bol(harvest_sendline,'+FIND_WIDTH'//NUL,tglf_find_width_in)
    ierr=set_harvest_payload_int(harvest_sendline,'+NWIDTH'//NUL,tglf_nwidth_in)

!   '#---------------------------------------------------'
!   '# Control parameters:'
!   '#---------------------------------------------------'
    ierr=set_harvest_payload_int(harvest_sendline,'+NS'//NUL,tglf_ns_in)
    ierr=set_harvest_payload_bol(harvest_sendline,'+ADIABATIC_ELEC'//NUL,tglf_adiabatic_elec_in)
    ierr=set_harvest_payload_bol(harvest_sendline,'+USE_BPER'//NUL,tglf_use_bper_in)
    ierr=set_harvest_payload_bol(harvest_sendline,'+USE_BPAR'//NUL,tglf_use_bpar_in)
    ierr=set_harvest_payload_bol(harvest_sendline,'+USE_MHD_RULE'//NUL,tglf_use_mhd_rule_in)
    ierr=set_harvest_payload_bol(harvest_sendline,'+USE_BISECTION'//NUL,tglf_use_bisection_in)
    ierr=set_harvest_payload_int(harvest_sendline,'+SAT_RULE'//NUL,tglf_sat_rule_in)
    ierr=set_harvest_payload_int(harvest_sendline,'+KYGRID_MODEL'//NUL,tglf_kygrid_model_in)
    ierr=set_harvest_payload_int(harvest_sendline,'+XNU_MODEL'//NUL,tglf_xnu_model_in)
    ierr=set_harvest_payload_int(harvest_sendline,'+VPAR_MODEL'//NUL,tglf_vpar_model_in)
    ierr=set_harvest_payload_int(harvest_sendline,'+VPAR_SHEAR_MODEL'//NUL,tglf_vpar_shear_model_in)
    ierr=set_harvest_payload_int(harvest_sendline,'+IBRANCH'//NUL,tglf_ibranch_in)
    ierr=set_harvest_payload_int(harvest_sendline,'+NMODES'//NUL,tglf_nmodes_in)
    ierr=set_harvest_payload_int(harvest_sendline,'+NXGRID'//NUL,tglf_nxgrid_in)
    ierr=set_harvest_payload_int(harvest_sendline,'+NKY'//NUL,tglf_nky_in)
    ierr=set_harvest_payload_dbl(harvest_sendline,'+ALPHA_MACH'//NUL,tglf_alpha_mach_in)
    ierr=set_harvest_payload_dbl(harvest_sendline,'+ALPHA_E'//NUL,tglf_alpha_e_in)
    ierr=set_harvest_payload_dbl(harvest_sendline,'+ALPHA_P'//NUL,tglf_alpha_p_in)
    ierr=set_harvest_payload_dbl(harvest_sendline,'+ALPHA_QUENCH'//NUL,tglf_alpha_quench_in)
    ierr=set_harvest_payload_dbl(harvest_sendline,'+XNU_FACTOR'//NUL,tglf_xnu_factor_in)
    ierr=set_harvest_payload_dbl(harvest_sendline,'+DEBYE_FACTOR'//NUL,tglf_debye_factor_in)
    ierr=set_harvest_payload_dbl(harvest_sendline,'+ETG_FACTOR'//NUL,tglf_etg_factor_in)

!   '#---------------------------------------------------'
!   '# Expert parameters:'
!   '#---------------------------------------------------'
    ierr=set_harvest_payload_dbl(harvest_sendline,'+DAMP_PSI'//NUL,tglf_damp_psi_in)
    ierr=set_harvest_payload_dbl(harvest_sendline,'+DAMP_SIG'//NUL,tglf_damp_sig_in)
    ierr=set_harvest_payload_dbl(harvest_sendline,'+PARK'//NUL,tglf_park_in)
    ierr=set_harvest_payload_dbl(harvest_sendline,'+GHAT'//NUL,tglf_ghat_in)
    ierr=set_harvest_payload_dbl(harvest_sendline,'+GCHAT'//NUL,tglf_gchat_in)
    ierr=set_harvest_payload_dbl(harvest_sendline,'+WD_ZERO'//NUL,tglf_wd_zero_in)
    ierr=set_harvest_payload_dbl(harvest_sendline,'+LINSKER_FACTOR'//NUL,tglf_linsker_factor_in)
    ierr=set_harvest_payload_dbl(harvest_sendline,'+GRADB_FACTOR'//NUL,tglf_gradB_factor_in)
    ierr=set_harvest_payload_dbl(harvest_sendline,'+FILTER'//NUL,tglf_filter_in)
    ierr=set_harvest_payload_dbl(harvest_sendline,'+THETA_TRAPPED'//NUL,tglf_theta_trapped_in)

!   '#---------------------------------------------------'
!   '# Disabled harvest:'
!   '#---------------------------------------------------'
!   ierr=set_harvest_payload_bol(harvest_sendline,'+USE_TRANSPORT_MODEL'//NUL,tglf_use_transport_model_in) ! DISABLED because we only allow runs with flux calculations
!   ierr=set_harvest_payload_int(harvest_sendline,'+GEOMETRY_FLAG'//NUL,tglf_geometry_flag_in) ! DISABLED because we only allow Miller geometry
!   ierr=set_harvest_payload_bol(harvest_sendline,'+NEW_EIKONAL'//NUL,tglf_new_eikonal_in) ! DISABLED because we only allow freshly computed eikonal
!   ierr=set_harvest_payload_dbl(harvest_sendline,'+VEXB'//NUL,tglf_vexb_in) !DISABLED because not in use, see VPAR
!   ierr=set_harvest_payload_bol(harvest_sendline,'+IFLUX'//NUL,tglf_iflux_in)  ! DISABLED because if not True then we do not come here
!   ierr=set_harvest_payload_int(harvest_sendline,'+NBASIS_MAX'//NUL,tglf_nbasis_max_in) ! DISABLED because not physics inputs/outputs
!   ierr=set_harvest_payload_int(harvest_sendline,'+NBASIS_MIN'//NUL,tglf_nbasis_min_in) ! DISABLED because not physics inputs/outputs
!   ierr=set_harvest_payload_int(harvest_sendline,'+WRITE_WAVEFUNCTION_FLAG'//NUL,tglf_write_wavefunction_flag_in)  ! DISABLED because not physics inputs/outputs

!   '#---------------------------------------------------'
!   '# Output results:'
!   '#---------------------------------------------------'
    DO i = 1,tglf_ns_in
      IF (i < 10) THEN
         write (NUM, "(I01,A1)") i,NUL
      else
         write (NUM, "(I02,A1)") i,NUL
      ENDIF
      ierr=set_harvest_payload_dbl(harvest_sendline,'OUT_PARTICLE_FLUX_'//NUM, &
                                   get_particle_flux(i,1)+get_particle_flux(i,2)+get_particle_flux(i,3))
      ierr=set_harvest_payload_dbl(harvest_sendline,'OUT_ENERGY_FLUX_'//NUM, &
                                   get_energy_flux(i,1)+get_energy_flux(i,2)+get_energy_flux(i,3))
      ierr=set_harvest_payload_dbl(harvest_sendline,'OUT_STRESS_PAR_'//NUM, &
                                   get_stress_par(i,1)+get_stress_par(i,2)+get_stress_par(i,3))
      ierr=set_harvest_payload_dbl(harvest_sendline,'OUT_STRESS_TOR_'//NUM, &
                                   get_stress_tor(i,1)+get_stress_tor(i,2)+get_stress_tor(i,3))
      ierr=set_harvest_payload_dbl(harvest_sendline,'OUT_STRESS_EXCHANGE_'//NUM, &
                                   get_exchange(i,1)+get_exchange(i,2)+get_exchange(i,3))
   ENDDO

   ierr=harvest_send(harvest_sendline)

  END SUBROUTINE tglf_harvest_local