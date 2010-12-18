!----------------------------------------------------------------
! neo_experimental_profiles.f90
!
! PURPOSE:
!  Read experimental profiles. and generate additional profiles.
!
! NOTES: 
!  See http://fusion.gat.com/theory/INPUT_profiles for a 
!  complete description of the INPUT_profiles file structure.
! 
!  The files read by this routine are *generated* from 
!  INPUT_profiles by scripts. 
!----------------------------------------------------------------

subroutine neo_experimental_profiles

  use neo_globals
  use neo_profile_exp
  use neo_allocate_profile
  use EXPRO_interface

  !---------------------------------------------------
  implicit none

  integer :: i, is, i_ion, i_exp

  !--------------------------------------------------------------
  ! use EXPRO routines to read data:
  !
  call EXPRO_alloc('./',1)
  EXPRO_ctrl_density_method = 2  ! quasi-neutrality density flag
  EXPRO_ctrl_signb = sign_bunit
  EXPRO_ctrl_signq = sign_q
  EXPRO_ctrl_rotation_method = 1 ! Standard method (Waltz=2 is not consistent with NEO).

  if (profile_equilibrium_model == 2) then
     EXPRO_ctrl_numeq_flag = 1
  else
     EXPRO_ctrl_numeq_flag = 0
  endif

  ! Determine if electrons are to be included in the neo 
  ! simulation.  Electron profiles are read even if not 
  ! to be included in the simulation (needed to re-scale 
  ! ion density/temp is not quasi-neutral).

  if (adiabatic_ele_model == 1) then
     n_species_exp = n_species + 1
  else
     if (Z(n_species) /= -1) then
        print *,&
             'For exp. profiles, electron species must be n_species'
        stop
     endif
     n_species_exp = n_species
  endif

  ! charge of ion species
  EXPRO_ctrl_z(:) = 0.0
  if(adiabatic_ele_model == 1) then
     do is=1,n_species
        EXPRO_ctrl_z(is) = 1.0 * Z(is)
     enddo
  else
     do is=1,n_species-1
        EXPRO_ctrl_z(is) = 1.0 * Z(is)
     enddo
  endif

  call EXPRO_read('./')

  n_grid_exp = EXPRO_n_exp
  call PROFILE_EXP_alloc(1)

  !--------------------------------------------------------------

  !--------------------------------------------------------------
  ! Transfer data from read vector to individual arrays:
  !
  rhoN_torflux_a      = EXPRO_arho
  rhoN_torflux_exp(:) = EXPRO_rho(:)
  psiN_polflux_exp(:) = EXPRO_poloidalfluxover2pi(:) / EXPRO_poloidalfluxover2pi(n_grid_exp)
  rmin_exp(:)         = EXPRO_rmin(:)
  rmaj_exp(:)         = EXPRO_rmaj(:)
  q_exp(:)            = EXPRO_q(:)
  kappa_exp(:)        = EXPRO_kappa(:)
  delta_exp(:)        = EXPRO_delta(:)  * profile_delta_scale
  zeta_exp(:)         = EXPRO_zeta(:)   * profile_zeta_scale
  zmag_exp(:)         = EXPRO_zmag(:)   * profile_zmag_scale
  shift_p(:)          = EXPRO_drmaj(:)
  s_kappa_p(:)        = EXPRO_skappa(:)
  s_delta_p(:)        = EXPRO_sdelta(:) * profile_delta_scale
  s_zeta_p(:)         = EXPRO_szeta(:)  * profile_zeta_scale
  s_zmag_p(:)         = EXPRO_dzmag(:)  * profile_zmag_scale
  shat_p(:)           = EXPRO_s(:)
  
  r_p(:)    = rmin_exp(:)/rmin_exp(n_grid_exp)
  rmaj_p(:) = rmaj_exp(:)/rmin_exp(n_grid_exp)
  zmag_p(:) = zmag_exp(:)/rmin_exp(n_grid_exp)

  ! Minor radius, a, in meters:
  a_meters = rmin_exp(n_grid_exp)

  if(write_out_mode > 1 .and. abs(profile_delta_scale-1.0) > epsilon(0.) ) then
     print *, 'Delta and S_Delta are re-scaled'
  end if
  if(write_out_mode > 1 .and. abs(profile_zeta_scale-1.0) > epsilon(0.)) then
     print *, 'Zeta and S_Zeta are re-scaled'
  end if
  if(write_out_mode > 1 .and. abs(profile_zmag_scale-1.0) > epsilon(0.)) then
     print *, 'Zmag and S_Zmag are re-scaled'
  end if

  if (profile_equilibrium_model == 2) then
     if(EXPRO_nfourier <= 0) then
        print *, 'ERROR: Geometry coefficients missing.'
        stop
     endif
     geo_numeq_flag = 1
     deallocate(geo_yin_exp)
     deallocate(geo_yin)
     geo_ny = EXPRO_nfourier
     allocate(geo_yin_exp(8,0:geo_ny,n_grid_exp))
     allocate(geo_yin(8,0:geo_ny,n_radial))
     geo_yin_exp(1:4,:,:) = EXPRO_geo(:,:,:)/rmin_exp(n_grid_exp)
     geo_yin_exp(5:8,:,:) = EXPRO_dgeo(:,:,:)
  endif

  te_ade_exp(:)     = EXPRO_te(:)
  ne_ade_exp(:)     = EXPRO_ne(:)

  tem_exp(n_species_exp,:)  = EXPRO_te(:)
  dlntdr_p(n_species_exp,:) = EXPRO_dlntedr(:) * a_meters
  den_exp(n_species_exp,:)  = EXPRO_ne(:)
  dlnndr_p(n_species_exp,:) = EXPRO_dlnnedr(:) * a_meters

  do i_ion=1,n_species_exp-1
     ! ion temps must be equal 
     if(profile_temprescale_model == 1) then
        if(write_out_mode > 1) then
           print *, 'Re-scaling ion temperatures for equal temps Ti=Te'
        end if
        tem_exp(i_ion,:)  = EXPRO_te(:)
        dlntdr_p(i_ion,:) = EXPRO_dlntedr(:) * a_meters
     else
        tem_exp(i_ion,:)  = EXPRO_ti(1,:)
        dlntdr_p(i_ion,:) = EXPRO_dlntidr(1,:) * a_meters
     endif
     ! first species density is re-set by quasi-neutrality
     if(i_ion == 1) then
        den_exp(i_ion,:)  = EXPRO_ni_new(:)
        dlnndr_p(i_ion,:) = EXPRO_dlnnidr_new(:) * a_meters
     else
        den_exp(i_ion,:)  = EXPRO_ni(i_ion,:)
        dlnndr_p(i_ion,:) = EXPRO_dlnnidr(i_ion,:) * a_meters
     endif
  enddo

  ! Sanity check for temperatures
  do i=1,n_species_exp
     if (minval(den_exp(i,:)) <= 0.0) then
        print *, 'ERROR: Nonpositive in exp. density profile'
        stop
     endif
  enddo

  ! Bunit 

  b_unit_p(:)  = EXPRO_bunit(:)

  ! Rotation and Er

  omega_rot_p(:)       = EXPRO_w0(:) 
  omega_rot_deriv_p(:) = EXPRO_w0p(:) 
  dphi0dr_p(:) = -omega_rot_p(:) * rmin_exp(:) * b_unit_p(:) / q_exp(:)

  call EXPRO_alloc('./',0)

end subroutine neo_experimental_profiles
