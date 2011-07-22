!----------------------------------------------------------------
! gyro_read_experimental_profiles.f90 [caller gyro_profile_init]
!
! PURPOSE:
!  Read experimental profiles using EXPRO library.
!
! NOTES: 
!  See http://fusion.gat.com/theory/input.profiles for a 
!  complete description of the input.profiles file structure.
!----------------------------------------------------------------

subroutine gyro_read_experimental_profiles

  use gyro_globals
  use gyro_profile_exp
  use EXPRO_interface

  !---------------------------------------------------
  implicit none
  !
  real :: p_total
  !---------------------------------------------------

  !---------------------------------------------------------------------
  ! Read experimental profiles using EXPRO library.
  !
  EXPRO_ctrl_density_method = density_method 
  EXPRO_ctrl_z = z_vec(1:5)
  EXPRO_ctrl_numeq_flag = num_equil_flag 
  EXPRO_ctrl_signq = ipccw*btccw
  EXPRO_ctrl_signb = -btccw
  EXPRO_ctrl_rotation_method = rotation_theory_method
  !
  ! Kludge 
  if (n_vec(2) == 0.0) EXPRO_ctrl_z(2) = 0.0
  if (n_vec(3) == 0.0) EXPRO_ctrl_z(3) = 0.0
  if (n_vec(4) == 0.0) EXPRO_ctrl_z(4) = 0.0
  if (n_vec(5) == 0.0) EXPRO_ctrl_z(5) = 0.0

  call EXPRO_palloc(GYRO_COMM_WORLD,path,1) 
  call EXPRO_pread

  if (i_proc == 0) call EXPRO_write_derived
  !---------------------------------------------------------------------

  !---------------------------------------------------------------------
  n_grid_exp = EXPRO_n_exp
  if (EXPRO_nfourier < 0) then
     n_fourier_geo = 0
     if (num_equil_flag == 1) then
        call catch_error('ERROR: Geometry coefficients missing.')
     endif
  else
     n_fourier_geo = EXPRO_nfourier
  endif
  !
  call allocate_profile_exp
  !
  ! Transfer data from read arrays to individual arrays:
  !
  bt_exp   = EXPRO_b_ref
  arho_exp = EXPRO_arho
  !
  rhogrid_exp(:)    = EXPRO_rho(:)
  rmin_exp(:)       = EXPRO_rmin(:)
  rmaj_exp(:)       = EXPRO_rmaj(:)
  zmag_exp(:)       = EXPRO_zmag(:)
  ptot_exp(:)       = EXPRO_ptot(:)
  q_exp(:)          = EXPRO_q(:)
  kappa_exp(:)      = EXPRO_kappa(:)
  delta_exp(:)      = EXPRO_delta(:)
  zeta_exp(:)       = EXPRO_zeta(:)
  tem_exp(n_spec,:) = EXPRO_te(:)
  den_exp(n_spec,:) = EXPRO_ne(:)
  z_eff_exp(:)      = EXPRO_z_eff(:)
  w0_exp(:)         = EXPRO_w0(:)
  w0p_exp(:)        = EXPRO_w0p(:)
  gamma_e_exp(:)    = EXPRO_gamma_e(:)
  gamma_p_exp(:)    = EXPRO_gamma_p(:)
  mach_exp(:)       = EXPRO_mach(:)
  !
  do is=1,n_spec-1
     den_exp(is,:)  = EXPRO_ni(is,:)*n_vec(is)
     tem_exp(is,:)  = EXPRO_ti(is,:)
     if (n_vec(is) /= 1.0) then
        call send_message_real(&
             'INFO: ni'//achar(is-1+iachar("1"))//' rescaled by: ',&
             n_vec(is))
     endif
  enddo

  !--------------------------------------------------------
  ! Minor radius, a, in meters:
  !
  a_meters = rmin_exp(n_grid_exp)
  !--------------------------------------------------------

  ! Fill in computed profile quantities:

  b_unit_p(:)    = EXPRO_bunit(:)
  shat_p(:)      = EXPRO_s(:)
  s_kappa_p(:)   = EXPRO_skappa(:)
  drmaj_p(:)     = EXPRO_drmaj(:)
  s_delta_p(:)   = EXPRO_sdelta(:)
  s_zeta_p(:)    = EXPRO_szeta(:)
  dzmag_p(:)     = EXPRO_dzmag(:)

  ! The case is=1 may be reset below.

  do is=1,n_spec-1
     dlntdr_p(is,:) = a_meters*EXPRO_dlntidr(is,:)
     dlnndr_p(is,:) = a_meters*EXPRO_dlnnidr(is,:)
  enddo
  dlntdr_p(n_spec,:) = a_meters*EXPRO_dlntedr(:)
  dlnndr_p(n_spec,:) = a_meters*EXPRO_dlnnedr(:)
  dlnptotdr_p(:)     = a_meters*EXPRO_dlnptotdr(:)
  !
  dlnndr_p(:,1)  = 0.0
  dlntdr_p(:,1)  = 0.0
  dlnptotdr_p(1) = 0.0
  !----------------------------------------------------------------------------

  ! Retain up-down asymmetry from elevation in the case of Miller shape
  if (num_equil_flag == 0 .and. udsymmetry_flag == 1) then
     zmag_exp(:) = 0.0
     dzmag_p(:)  = 0.0
  endif

  ! Fill in general geometry parameters if they exist
  if (num_equil_flag == 1) then

     ! Fill in geometry arrays with original data from 
     ! EXPRO routine.  Note that EXPRO_geo is in m, 
     ! and EXPRO_dgeo is dimensionless.  All geo_p are 
     ! dimensionless:

     geo_p(1:4,:,:) = EXPRO_geo(:,:,:)/a_meters
     geo_p(5:8,:,:) = EXPRO_dgeo(:,:,:)

  endif

  ! rhos_deuterium/a (dimensionless)
  rhosda_p(:) = EXPRO_rhos(:)/a_meters

  ! cs_deuterium/a (1/s)
  csda_p(:)   = EXPRO_cs(:)/a_meters

  !--------------------------------------------------------------
  ! Management of multiple ion densities:
  !
  select case (density_method)

  case (1)

     call send_message('INFO: Taking densities directly from INPUT_profiles')

  case (2) 

     call send_message('INFO: Offsetting main ions to force sum(ni) = ne')
     den_exp(1,:) = EXPRO_ni_new(:)
     dlnndr_p(1,:) = a_meters*EXPRO_dlnnidr_new(:)

  case default

     call catch_error('INVALID: density_method')

  end select
  !--------------------------------------------------------------

  !-------------------------------------------------------------------
  ! Sanity check for densities
  !
  do is=1,n_spec
     if (minval(den_exp(is,:)) <= 0.0) then
        print *,is,den_exp(is,:)
        call catch_error(&
             'ERROR: Nonpositive in density profile '//achar(is-1+iachar("1")))  
     endif
  enddo

  do i_exp=1,n_grid_exp

     r_p(i_exp) = rmin_exp(i_exp)/rmin_exp(n_grid_exp)
     q_p(i_exp) = q_exp(i_exp)

     if (q_scale /= 1.0) then
        q_p(i_exp) = q_exp(i_exp)*q_scale
     endif

     !-------------------------------------------------------------
     ! Specification of beta_unit (beta in terms of b_unit).
     !
     ! beta_unit is defined with respect to b_unit, as required 
     ! by GEO.  For comparison with DIII-D, beta should be 
     ! defined with respect to bt_exp.
     !
     p_total = sum(den_exp(:,i_exp)*tem_exp(:,i_exp))

     beta_unit_p(i_exp) = 400.0*p_total/(1e5*b_unit_p(i_exp)**2)
     beta_unit_ptot_p(i_exp) = 400.0*ptot_exp(i_exp)/(1.6022*1e3)/(1.e5*b_unit_p(i_exp)**2)
     !-------------------------------------------------------------

  enddo

  if (ampere_scale /= 1.0) then
     call send_message_real('INFO: betae in Ampere Eq. scaled by ',ampere_scale)
  endif

  ! Deallocate EXPRO arrays
  call EXPRO_palloc(GYRO_COMM_WORLD,path,0)

  if (debug_flag == 1 .and. i_proc == 0) then
     print *,'[read_experimental_profiles done]'
  endif
  !---------------------------------------------------------------------

end subroutine gyro_read_experimental_profiles
