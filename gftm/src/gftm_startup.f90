       SUBROUTINE gftm_startup
!*********************************************
!
!     initalization of gftm
!
!*********************************************
      USE gftm_dimensions
      USE gftm_global
      USE gftm_velocity_matrix
      USE gftm_gyro_average
!
      IMPLICIT NONE
      INTEGER :: i,j
      LOGICAL :: USE_PRESETS = .TRUE.
!
! initialized trace_path that records the flow of gftm internal calls
!
      do i=1,7
        trace_path(i)=0
      enddo
!
! set global constants
!
      pi = atan2( 0.0, -1.0 )
      pi_2 = 2.0*pi
      sqrt_pi = SQRT(pi)
      sqrt_two = SQRT(2.0)
!
! set dimensions of allocatable arrays
!
      ns0 = 1
      if(adiabatic_elec_in) then
        ns0 = 2
      endif
      if(use_default_species)ns_in=2
      ns = ns_in
!      write(*,*)"ns = ",ns,"   ns0 = ",ns0
      nx = 2*nxgrid_in-1
      nbasis_max_in = 8
      nbasis = nbasis_max_in
      nbasis_max = nbasis
! case 1
!      nu = 3
!      ne = 2
!      damp_psi_in = 0.06
!      damp_sig_in = 0.15
! case 2
      nu = 5
      ne = 3
 !     damp_psi_in = 0.028
 !     damp_sig_in = 0.18
 !     damp_psi_in = 0.015
 !     damp_sig_in = 0.15
 ! for width = 1.65, nbasis = 6
 !     damp_psi_in = 0.0075
 !     damp_sig_in = 0.083
 ! for width = 1.75, nbasis = 6
 !     damp_psi_in = 0.020
 !     damp_sig_in = 0.15
 !     etg_factor_in = 0.391827
 !     damp_psi_in = 0.0
 !     damp_sig_in = 0.5
 !     etg_factor_in = 0.33
 !     damp_psi_in = 0.0
 !      damp_sig_in = 0.75
 !      etg_factor_in = 0.25
!
! case 3
 !     nu = 7
 !     ne = 4
!      damp_psi_in = 0.006
 !     damp_sig_in = 0.12
! case 4
!      nu = 12
!      ne = 6
!      damp_psi_in = 0.005
!      damp_sig_in = 0.09
!
!   values used for 5-9
!      damp_psi_in = 0.005
!      damp_sig_in = 0.087
!
      nune = nu*ne
      nphase = nu*ne*nbasis
      ntot = (ns-ns0+1)*nu*ne*nbasis
!
      ! restrict user settings to three preset versions
      ! SAT0 + GYRO + XNU 2
      ! SAT1 + XNU 2
      ! SAT2 + XNU 3 + UNITS != GYRO
      if(USE_PRESETS)then
        wdia_trapped_in=0.0
        if(sat_rule_in.eq.2)then
          xnu_model_in=3
          wdia_trapped_in = 1.0
          if(units_in.eq."GYRO")units_in = "CGYRO"
          if(igeo .ne. 1)then
            write(*,*)"SAT_RULE=2 requires Miller geometry, set GEOMETRY_FLAG=1"
            STOP
          endif
        endif
        if(sat_rule_in.eq.1)then
          xnu_model_in = 2
        endif
        if(sat_rule_in.eq.0)then
          units_in = 'GYRO'
          xnu_model_in = 2
        endif
        if(use_bper_in)alpha_mach_in=0.0
      endif
! SAT_RULE=0 has different fits for nmodes=2 and nmodes=4
      if(sat_rule_in.eq.0)then
        if(nmodes_in.gt.2)nmodes_in=4
      endif
! debug
!      write(*,*)"nx=",nx
!      write(*,*)"nbasis =",nbasis    
!      write(*,*)"ns0 = ",ns0
!      write(*,*)"ns = ",ns
!      write(*,*)"nxgrid_in=",nxgrid_in
!      write(*,*)"nky_in=",nky_in
!      write(*,*)"nmodes_in=",nmodes_in
!      write(*,*)"nbasis_max_in=",nbasis_max_in
!      write(*,*)"nbasis_min_in=",nbasis_min_in
!      write(*,*)"nwidth_in=",nwidth_in
!
! allocate internal matricies
!
      CALL gftm_allocate
!
! fill species arrays
!
        trace_path(3)=1
      CALL get_species
!      write(*,*)"nxgrid_in=",nxgrid_in
!      write(*,*)"nky_in=",nky_in
!      write(*,*)"nmodes_in=",nmodes_in
!      write(*,*)"nbasis_max_in=",nbasis_max_in
!      write(*,*)"nbasis_min_in=",nbasis_min_in
!      write(*,*)"nwidth_in=",nwidth_in
!
      if(new_start)then
        trace_path(2)=1
        new_width=.TRUE.
      endif
      new_start = .FALSE.
!
!   load the xgrid and gauss-hermite weights
!
      if(gauher_uncalled)CALL gauher
!
! load the hermite basis functions on the xgrid
!  
      if(gauss_hermite_uncalled)CALL gauss_hermite
!
      if(igeo.eq.1)grad_r0_out = 1.0/(1.0 + drmajdx_loc)
!
!debug
      go to 10
      call get_velocity_matrix
      write(*,*)"mat_upar"
      do i=1,num
       write(*,*)(mat_upar(i,j),j=1,num)
      enddo
      write(*,*)"mat_dupar"
      do i=1,num
       write(*,*)(mat_dupar(i,j),j=1,num)
      enddo
      write(*,*)"mat_eper"
      do i=1,nem
       write(*,*)(mat_eper(i,j),j=1,nem)
      enddo
      write(*,*)"mat_deper"
      do i=1,nem
       write(*,*)(mat_deper(i,j),j=1,nem)
      enddo
      call get_gyro_average(1.0,nem)
      write(*,*)"pe1j0",(pe1j0(j),j=1,nem)
      write(*,*)"pe2j0",(pe2j0(j),j=1,nem)
      write(*,*)"pe1j1",(pe1j1(j),j=1,nem)
      write(*,*)"pe2j1",(pe2j1(j),j=1,nem)
 10   continue
! debug
      END SUBROUTINE gftm_startup
!
!
      SUBROUTINE get_species
!******************************************************************************
!  assign species arrays
!******************************************************************************
      USE gftm_dimensions
      USE gftm_global
      USE gftm_species
!
      IMPLICIT NONE
      INTEGER :: is
      REAL :: charge
!
! electrons=1, ions =2,...
!
      if(use_default_species)then
       mass_in(1) = 0.0002723
       mass_in(2) = 1.0
       zs_in(1) = -1.0
       zs_in(2) = 1.0
       as_in(1) = 1.0
       as_in(2) = 1.0
       taus_in(1)=1.0
       taus_in(2)=1.0
       ns_in = 2
       nstotal_in = 2
      endif
!
      nfields_out = 1
      if(use_bper_in)nfields_out = nfields_out + 1
      if(use_bpar_in)nfields_out = nfields_out + 1
!  inputs
      ky_s = ky_in
      vexb_shear_s = vexb_shear_in*sign_It_in
      xnue_s = xnue_in
      pol = 0.0
      U0 = 0.0
      charge = 0.0
      rho_ion = 0.0
!      do is=1,nstotal_in  ! include all species inputs
      do is=1,ns
        rlns(is) = rlns_in(is)
        rlts(is) = rlts_in(is)
        taus(is) = taus_in(is)
        as(is) = as_in(is)
        zs(is) = zs_in(is)
        mass(is) = mass_in(is)
        vpar_s(is) = alpha_mach_in*sign_It_in*vpar_in(is)
        vpar_shear_s(is) = alpha_p_in*sign_It_in*vpar_shear_in(is)
        vs(is) = SQRT(taus(is)/mass(is))
        pol = pol +  zs(is)*zs(is)*as(is)/taus(is)
        U0 = U0 + as(is)*vpar_s(is)*zs(is)*zs(is)/taus(is)
        fts(is) = 0.0
        if(is.gt.1.and.zs(is)*as(is)/ABS(as(1)*zs(1)).gt.0.1)then
          charge = charge + zs(is)*as(is)
          rho_ion = rho_ion + zs(is)*as(is)*SQRT(2.0*mass(is)*taus(is))/zs(is) ! charge weighted average ion gyroradius
        endif
        rho_e =SQRT(2.0*mass(1)*taus(1))/ABS(zs(1))
!        write(*,*)"species",is
!        write(*,*)" vs = ",vs(is)
!        write(*,*)rlns(is),rlts(is)
!        write(*,*)"taus = ",taus(is),"   mass = ",mass(is)
!        write(*,*)zs(is),as(is)
      enddo
      if(charge.lt.0.001)call gftm_error(1,"total ion charge < 0.001")
      rho_ion = rho_ion/charge
      if(use_ave_ion_grid_in .eqv. .false.)then
        rho_ion = SQRT(2.0*mass(2)*taus(2))/zs(2)
      endif
!      write(*,*)"rho_ion = ",rho_ion
!      write(*,*)"rho_e = ",rho_e
!      write(*,*)"charge = ",charge
      END SUBROUTINE get_species
