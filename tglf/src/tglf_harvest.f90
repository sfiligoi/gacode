  ! Harvest LOCAL INTERFACE variables
  !
  ! _bol : boolean
  ! _int : integer
  ! _dbl : double
  !
  ! NOTE: here we refer to quantities defined in the write_tglf_input subroutine of tglf_inout.f90
  !       because that is actually what is used by TGLF at runtime

  SUBROUTINE tglf_harvest_local()

    USE tglf_interface
    USE tglf_pkg
    USE tglf_global

    INTEGER :: ierr, i, j, nky
    CHARACTER(LEN=65507) :: harvest_sendline
    CHARACTER(LEN=255) :: harvest_tag
    CHARACTER(LEN=2) :: NUM
    CHARACTER NUL
    PARAMETER(NUL = CHAR(0))

    REAL, DIMENSION(10) :: tmp
    INTEGER, DIMENSION(10) :: ions_order
    
    REAL, DIMENSION(:), ALLOCATABLE :: spectrum

    IF (.NOT.tglf_use_transport_model_in) THEN
        WRITE(1,*) 'HARVEST ONLY WHEN `TRANSPORT_MODEL=.TRUE.`' !only when computing fluxes
        RETURN
    ENDIF

    IF (.NOT.iflux_in) THEN
        WRITE(1,*) 'HARVEST ONLY AVAILABLE WHEN `IFLUX=.TRUE.' !only when computing fluxes
        RETURN
    ENDIF

    IF (igeo .NE. 1 ) THEN
       WRITE(1,*) 'HARVEST ONLY SUPPORTS MILLER `GEOMETRY_FLAG=1`' !only when using miller
       RETURN
    ENDIF

    ierr=init_harvest('TGLF_spectrum_3'//NUL,harvest_sendline,LEN(harvest_sendline))

    ierr=set_harvest_payload_str(harvest_sendline,'VERSION'//NUL,'APS15_1'//NUL) !no underscore to allow different versions of the same run

!   ---------------------------------------------------
!    Plasma parameters
!   ---------------------------------------------------
    ierr=set_harvest_payload_int(harvest_sendline,'SIGN_BT'//NUL,sign_Bt_in)
    ierr=set_harvest_payload_int(harvest_sendline,'SIGN_IT'//NUL,sign_It_in)
    IF (kygrid_model_in.NE.1) THEN
        ierr=set_harvest_payload_dbl(harvest_sendline,'KY'//NUL,ky_in)
    ENDIF
    ierr=set_harvest_payload_dbl(harvest_sendline,'VEXB_SHEAR'//NUL,vexb_shear_in)
    ierr=set_harvest_payload_dbl(harvest_sendline,'BETAE'//NUL,betae_in)
    ierr=set_harvest_payload_dbl(harvest_sendline,'XNUE'//NUL,xnue_in)
    ierr=set_harvest_payload_dbl(harvest_sendline,'ZEFF'//NUL,zeff_in)
    IF (debye_in.NE.0) THEN
        ierr=set_harvest_payload_dbl(harvest_sendline,'DEBYE'//NUL,debye_in)
    ENDIF

!   ---------------------------------------------------
!    Sort ions by A, Z, Te/Ti
!    electrons always first (as required by TGLF)
!    !main ion always second (affects GB normalization)
!   ---------------------------------------------------
    tmp=0.0
    tmp(1)=1E10
    tmp(2)=1E9
    ions_order=0
    DO i = 3,ns_in
        tmp(i)=mass_in(i)*100000+zs_in(i)*1000+1./(1+rlts_in(i))
    ENDDO
    DO i = 1,ns_in
        j=MAXLOC(tmp, DIM=1)
        ions_order(i)=j
        tmp(j)=0
    ENDDO

!   ---------------------------------------------------
!    Species vectors
!   ---------------------------------------------------
    DO i = 1,ns_in
      IF (i < 10) THEN
         write (NUM, "(I01,A1)") i,NUL
      ELSE
         write (NUM, "(I02,A1)") i,NUL
      ENDIF

      j=ions_order(i)

      ierr=set_harvest_payload_dbl(harvest_sendline,'+ZS_'//NUM,zs_in(j))     !make new tables for different charges
      ierr=set_harvest_payload_dbl(harvest_sendline,'+MASS_'//NUM,mass_in(j)) !make new tables for different masses
      ierr=set_harvest_payload_dbl(harvest_sendline,'RLNS_'//NUM,rlns_in(j))
      ierr=set_harvest_payload_dbl(harvest_sendline,'RLTS_'//NUM,rlts_in(j))
      if (i.NE.1) THEN !these are always 1, by definition
        ierr=set_harvest_payload_dbl(harvest_sendline,'TAUS_'//NUM,taus_in(j))
        ierr=set_harvest_payload_dbl(harvest_sendline,'AS_'//NUM,as_in(j))
      ENDIF
      ierr=set_harvest_payload_dbl(harvest_sendline,'VPAR_'//NUM,vpar_in(j))
      ierr=set_harvest_payload_dbl(harvest_sendline,'VPAR_SHEAR_'//NUM,vpar_shear_in(j))
    ENDDO

!   ---------------------------------------------------
!    Miller geometry parameters
!   ---------------------------------------------------
    ierr=set_harvest_payload_dbl(harvest_sendline,'RMIN_LOC'//NUL,rmin_loc)
    ierr=set_harvest_payload_dbl(harvest_sendline,'RMAJ_LOC'//NUL,rmaj_loc)
    ierr=set_harvest_payload_dbl(harvest_sendline,'ZMAJ_LOC'//NUL,zmaj_loc)
    IF (drmindx_loc.NE.1) THEN
        ierr=set_harvest_payload_dbl(harvest_sendline,'DRMINDX_LOC'//NUL,drmindx_loc) ! different derivatives affect everything!
    ENDIF
    ierr=set_harvest_payload_dbl(harvest_sendline,'DRMAJDX_LOC'//NUL,drmajdx_loc)
    ierr=set_harvest_payload_dbl(harvest_sendline,'DZMAJDX_LOC'//NUL,dzmajdx_loc)
    ierr=set_harvest_payload_dbl(harvest_sendline,'Q_LOC'//NUL,q_loc)
    ierr=set_harvest_payload_dbl(harvest_sendline,'KAPPA_LOC'//NUL,kappa_loc)
    ierr=set_harvest_payload_dbl(harvest_sendline,'S_KAPPA_LOC'//NUL,s_kappa_loc)
    ierr=set_harvest_payload_dbl(harvest_sendline,'DELTA_LOC'//NUL,delta_loc)
    ierr=set_harvest_payload_dbl(harvest_sendline,'S_DELTA_LOC'//NUL,s_delta_loc)
    ierr=set_harvest_payload_dbl(harvest_sendline,'ZETA_LOC'//NUL,zeta_loc)
    ierr=set_harvest_payload_dbl(harvest_sendline,'S_ZETA_LOC'//NUL,s_zeta_loc)
    ierr=set_harvest_payload_dbl(harvest_sendline,'P_PRIME_LOC'//NUL,p_prime_loc)
    ierr=set_harvest_payload_dbl(harvest_sendline,'Q_PRIME_LOC'//NUL,q_prime_loc)
    IF (kx0_loc.NE.0) THEN
        ierr=set_harvest_payload_dbl(harvest_sendline,'KX0_LOC'//NUL,kx0_loc) !should be always 0.0
    ENDIF

!   ---------------------------------------------------
!    Gaussian width parameters
!   ---------------------------------------------------
    ierr=set_harvest_payload_dbl(harvest_sendline,'+WIDTH'//NUL,width_in)
    ierr=set_harvest_payload_dbl(harvest_sendline,'+WIDTH_MIN'//NUL,width_min_in)
    ierr=set_harvest_payload_bol(harvest_sendline,'+FIND_WIDTH'//NUL,find_width_in)
    ierr=set_harvest_payload_int(harvest_sendline,'+NWIDTH'//NUL,nwidth_in)

!   ---------------------------------------------------
!    Control parameters
!   ---------------------------------------------------
    ierr=set_harvest_payload_bol(harvest_sendline,'+USE_TRANSPORT_MODEL'//NUL,tglf_use_transport_model_in)
    ierr=set_harvest_payload_int(harvest_sendline,'+GEOMETRY_FLAG'//NUL,igeo)
    ierr=set_harvest_payload_bol(harvest_sendline,'+IFLUX'//NUL,iflux_in)
    ierr=set_harvest_payload_int(harvest_sendline,'+NS'//NUL,ns_in)
    ierr=set_harvest_payload_bol(harvest_sendline,'+ADIABATIC_ELEC'//NUL,adiabatic_elec_in)
    ierr=set_harvest_payload_bol(harvest_sendline,'+USE_BPER'//NUL,use_bper_in)
    ierr=set_harvest_payload_bol(harvest_sendline,'+USE_BPAR'//NUL,use_bpar_in)
    ierr=set_harvest_payload_bol(harvest_sendline,'+USE_MHD_RULE'//NUL,use_mhd_rule_in)
    ierr=set_harvest_payload_bol(harvest_sendline,'+USE_BISECTION'//NUL,use_bisection_in)
    ierr=set_harvest_payload_int(harvest_sendline,'+SAT_RULE'//NUL,sat_rule_in)
    ierr=set_harvest_payload_int(harvest_sendline,'+KYGRID_MODEL'//NUL,kygrid_model_in)
    ierr=set_harvest_payload_int(harvest_sendline,'+XNU_MODEL'//NUL,xnu_model_in)
    ierr=set_harvest_payload_int(harvest_sendline,'+VPAR_MODEL'//NUL,vpar_model_in)
    ierr=set_harvest_payload_int(harvest_sendline,'+VPAR_SHEAR_MODEL'//NUL,vpar_shear_model_in)
    ierr=set_harvest_payload_int(harvest_sendline,'+IBRANCH'//NUL,ibranch_in)
    ierr=set_harvest_payload_int(harvest_sendline,'+NMODES'//NUL,nmodes_in)
    ierr=set_harvest_payload_int(harvest_sendline,'+NBASIS_MAX'//NUL,nbasis_max_in)
    ierr=set_harvest_payload_int(harvest_sendline,'+NBASIS_MIN'//NUL,nbasis_min_in)
    ierr=set_harvest_payload_int(harvest_sendline,'+NXGRID'//NUL,nxgrid_in)
    ierr=set_harvest_payload_int(harvest_sendline,'+NKY'//NUL,nky_in)
    ierr=set_harvest_payload_dbl(harvest_sendline,'+ALPHA_MACH'//NUL,alpha_mach_in)
    ierr=set_harvest_payload_dbl(harvest_sendline,'+ALPHA_E'//NUL,alpha_e_in)
    ierr=set_harvest_payload_dbl(harvest_sendline,'+ALPHA_P'//NUL,alpha_p_in)
    ierr=set_harvest_payload_dbl(harvest_sendline,'+ALPHA_QUENCH'//NUL,alpha_quench_in)
    ierr=set_harvest_payload_dbl(harvest_sendline,'+ALPHA_ZF'//NUL,alpha_zf_in)
    ierr=set_harvest_payload_dbl(harvest_sendline,'+XNU_FACTOR'//NUL,xnu_factor_in)
    ierr=set_harvest_payload_dbl(harvest_sendline,'+DEBYE_FACTOR'//NUL,debye_factor_in)
    ierr=set_harvest_payload_dbl(harvest_sendline,'+ETG_FACTOR'//NUL,etg_factor_in)

!   ---------------------------------------------------
!    Expert parameters
!   ---------------------------------------------------
    ierr=set_harvest_payload_dbl(harvest_sendline,'+DAMP_PSI'//NUL,damp_psi_in)
    ierr=set_harvest_payload_dbl(harvest_sendline,'+DAMP_SIG'//NUL,damp_sig_in)
    ierr=set_harvest_payload_dbl(harvest_sendline,'+PARK'//NUL,park_in)
    ierr=set_harvest_payload_dbl(harvest_sendline,'+GHAT'//NUL,ghat_in)
    ierr=set_harvest_payload_dbl(harvest_sendline,'+GCHAT'//NUL,gchat_in)
    ierr=set_harvest_payload_dbl(harvest_sendline,'+WD_ZERO'//NUL,wd_zero_in)
    ierr=set_harvest_payload_dbl(harvest_sendline,'+LINSKER_FACTOR'//NUL,Linsker_factor_in)
    ierr=set_harvest_payload_dbl(harvest_sendline,'+GRADB_FACTOR'//NUL,gradB_factor_in)
    ierr=set_harvest_payload_dbl(harvest_sendline,'+FILTER'//NUL,filter_in)
    ierr=set_harvest_payload_dbl(harvest_sendline,'+THETA_TRAPPED'//NUL,theta_trapped_in)

!   ---------------------------------------------------
!    Disabled harvest
!   ---------------------------------------------------
!   ierr=set_harvest_payload_bol(harvest_sendline,'+NEW_EIKONAL'//NUL,tglf_new_eikonal_in) !DISABLED because not used anymore and should not make difference anyways
!   ierr=set_harvest_payload_dbl(harvest_sendline,'+VEXB'//NUL,vexb_in) !DISABLED because not used anymore
!   ierr=set_harvest_payload_int(harvest_sendline,'+WRITE_WAVEFUNCTION_FLAG'//NUL,tglf_write_wavefunction_flag_in)  ! DISABLED because it does not matter for us
!   ierr=set_harvest_payload_dbl(harvest_sendline,'VNS_SHEAR_'//NUM,vns_shear_in(j)) !DISABLED because not used anymore
!   ierr=set_harvest_payload_dbl(harvest_sendline,'VTS_SHEAR_'//NUM,vts_shear_in(j)) !DISABLED because not used anymore

!   ---------------------------------------------------
!    Output results
!   ---------------------------------------------------
    DO i = 1,ns_in
      IF (i < 10) THEN
         write (NUM, "(I01,A1)") i,NUL
      ELSE
         write (NUM, "(I02,A1)") i,NUL
      ENDIF

      j=ions_order(i)

      ierr=set_harvest_payload_dbl(harvest_sendline,'OUT_PARTICLE_FLUX_ESTATIC_'//NUM,get_particle_flux(j,1))
      ierr=set_harvest_payload_dbl(harvest_sendline,'OUT_ENERGY_FLUX_ESTATIC_'//NUM,get_energy_flux(j,1))
      ierr=set_harvest_payload_dbl(harvest_sendline,'OUT_STRESS_PAR_ESTATIC_'//NUM,get_stress_par(j,1))
      ierr=set_harvest_payload_dbl(harvest_sendline,'OUT_STRESS_TOR_ESTATIC_'//NUM,get_stress_tor(j,1))
      ierr=set_harvest_payload_dbl(harvest_sendline,'OUT_STRESS_EXCHANGE_ESTATIC_'//NUM,get_exchange(j,1))

      if (use_bper_in) THEN
          ierr=set_harvest_payload_dbl(harvest_sendline,'OUT_PARTICLE_FLUX_EMPER_'//NUM,get_particle_flux(j,2))
          ierr=set_harvest_payload_dbl(harvest_sendline,'OUT_ENERGY_FLUX_EMPER_'//NUM,get_energy_flux(j,2))
          ierr=set_harvest_payload_dbl(harvest_sendline,'OUT_STRESS_PAR_EMPER_'//NUM,get_stress_par(j,2))
          ierr=set_harvest_payload_dbl(harvest_sendline,'OUT_STRESS_TOR_EMPER_'//NUM,get_stress_tor(j,2))
          ierr=set_harvest_payload_dbl(harvest_sendline,'OUT_STRESS_EXCHANGE_EMPER_'//NUM,get_exchange(j,2))
      ENDIF

      if (use_bpar_in) THEN
          ierr=set_harvest_payload_dbl(harvest_sendline,'OUT_PARTICLE_FLUX_EMPAR_'//NUM,get_particle_flux(j,3))
          ierr=set_harvest_payload_dbl(harvest_sendline,'OUT_ENERGY_FLUX_EMPAR_'//NUM,get_energy_flux(j,3))
          ierr=set_harvest_payload_dbl(harvest_sendline,'OUT_STRESS_PAR_EMPAR_'//NUM,get_stress_par(j,3))
          ierr=set_harvest_payload_dbl(harvest_sendline,'OUT_STRESS_TOR_EMPAR_'//NUM,get_stress_tor(j,3))
          ierr=set_harvest_payload_dbl(harvest_sendline,'OUT_STRESS_EXCHANGE_EMPAR_'//NUM,get_exchange(j,3))
      ENDIF

   ENDDO

!   ---------------------------------------------------
!    Spectra
!   ---------------------------------------------------
   nky = get_nky_out()

   ALLOCATE(spectrum(nky))
   
   DO i = 1, nky
      spectrum(i) = get_ky_spectrum_out(i)
   ENDDO

   ierr=set_harvest_payload_dbl_array(harvest_sendline,'KY_SPECTRUM'//NUL,spectrum,nky)
   
   DO i = 1, nmodes_in
      IF (i < 10) THEN
         write (NUM, "(I01,A1)") i,NUL
      ELSE
         write (NUM, "(I02,A1)") i,NUL
      ENDIF
      
      DO j = 1, nky
         spectrum(j) = get_eigenvalue_spectrum_out(1,j,i)
      ENDDO 
      ierr=set_harvest_payload_dbl_array(harvest_sendline,'OUT_EIGENVALUE_SPECTRUM_GAMMA'//NUM,spectrum,nky)
      DO j = 1, nky
         spectrum(j) = get_eigenvalue_spectrum_out(2,j,i)
      ENDDO
      ierr=set_harvest_payload_dbl_array(harvest_sendline,'OUT_EIGENVALUE_SPECTRUM_OMEGA'//NUM,spectrum,nky)
   ENDDO
   DEALLOCATE(spectrum)

!   ---------------------------------------------------
!    Additional entries only if an harvest_tag is set
!   ---------------------------------------------------
   harvest_tag = NUL
   ierr = get_harvest_tag(harvest_tag,len(harvest_tag))
   if ( len_trim(harvest_tag) .gt. 1) then
      ierr=set_harvest_payload_raw(harvest_sendline,TRIM(tglf_harvest_extra_in)//NUL)
   endif

!   ---------------------------------------------------
!    Send data
!   ---------------------------------------------------
   ierr=harvest_send(harvest_sendline)

  END SUBROUTINE tglf_harvest_local
