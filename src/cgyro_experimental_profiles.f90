!----------------------------------------------------------------
! cgyro_experimental_profiles.f90
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

subroutine cgyro_experimental_profiles

  use cgyro_globals
  use cgyro_profile_exp
  use cgyro_allocate_profile
  use EXPRO_interface

  !---------------------------------------------------
  implicit none

  integer :: i, is, i_ion, i_exp

  !--------------------------------------------------------------
  ! use EXPRO routines to read data:
  !
  call EXPRO_alloc(path,1)
  EXPRO_ctrl_density_method = 2  ! quasi-neutrality density flag
  EXPRO_ctrl_rotation_method = 1 ! Standard method 

  if (equilibrium_model == 3) then
     EXPRO_ctrl_numeq_flag = 1
  else
     EXPRO_ctrl_numeq_flag = 0
  endif

  ! Determine if electrons are to be included in the cgyro 
  ! simulation.  Electron profiles are read even if not 
  ! to be included in the simulation (needed to re-scale 
  ! ion density/temp is not quasi-neutral).

  if (adiabatic_ele_model == 1) then
     n_species_exp = n_species + 1
  else
     if (z(n_species) /= -1) then
        call cgyro_error('ERROR: (GKCOLL) For exp. profiles, electron species must be n_species')
        return
     endif
     n_species_exp = n_species
  endif

  ! charge of ion species
  EXPRO_ctrl_z(:) = 0.0
  if(adiabatic_ele_model == 1) then
     do is=1,n_species
        EXPRO_ctrl_z(is) = 1.0 * z(is)
     enddo
  else
     do is=1,n_species-1
        EXPRO_ctrl_z(is) = 1.0 * z(is)
     enddo
  endif

  call EXPRO_read

  n_grid_exp = EXPRO_n_exp
  call PROFILE_EXP_alloc(1)

  sign_bunit = EXPRO_signb
  sign_q     = EXPRO_signq

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
  delta_exp(:)        = EXPRO_delta(:)  
  zeta_exp(:)         = EXPRO_zeta(:)   
  zmag_exp(:)         = EXPRO_zmag(:)   
  shift_p(:)          = EXPRO_drmaj(:)
  s_kappa_p(:)        = EXPRO_skappa(:)
  s_delta_p(:)        = EXPRO_sdelta(:) 
  s_zeta_p(:)         = EXPRO_szeta(:)  
  s_zmag_p(:)         = EXPRO_dzmag(:) 
  shat_p(:)           = EXPRO_s(:)
  
  r_p(:)    = rmin_exp(:)/rmin_exp(n_grid_exp)
  rmaj_p(:) = rmaj_exp(:)/rmin_exp(n_grid_exp)
  zmag_p(:) = zmag_exp(:)/rmin_exp(n_grid_exp)

  ! Minor radius, a, in meters:
  a_norm = rmin_exp(n_grid_exp)

  if (equilibrium_model == 3) then
     if(EXPRO_nfourier <= 0) then
        call cgyro_error('ERROR: (GKCOLL) Geometry coefficients missing')
        return
     endif
     geo_numeq_flag = 1
     deallocate(geo_yin_exp)
     deallocate(geo_yin)
     geo_ny = EXPRO_nfourier
     allocate(geo_yin_exp(8,0:geo_ny,n_grid_exp))
     allocate(geo_yin(8,0:geo_ny))
     geo_yin_exp(1:4,:,:) = EXPRO_geo(:,:,:)/rmin_exp(n_grid_exp)
     geo_yin_exp(5:8,:,:) = EXPRO_dgeo(:,:,:)
  endif

  te_ade_exp(:)     = EXPRO_te(:)
  ne_ade_exp(:)     = EXPRO_ne(:)

  tem_exp(n_species_exp,:)  = EXPRO_te(:)
  dlntdr_p(n_species_exp,:) = EXPRO_dlntedr(:) * a_norm
  den_exp(n_species_exp,:)  = EXPRO_ne(:)
  dlnndr_p(n_species_exp,:) = EXPRO_dlnnedr(:) * a_norm

  do i_ion=1,n_species_exp-1
     tem_exp(i_ion,:)  = EXPRO_ti(i_ion,:)
     dlntdr_p(i_ion,:) = EXPRO_dlntidr(i_ion,:) * a_norm
     den_exp(i_ion,:)  = EXPRO_ni(i_ion,:)
     dlnndr_p(i_ion,:) = EXPRO_dlnnidr(i_ion,:) * a_norm
  enddo

  ! Sanity check for density
  do i=1,n_species_exp
     if (minval(den_exp(i,:)) <= 0.0) then
        call cgyro_error('ERROR: (GKCOLL) Nonpositive in exp. density profile')
        return
     endif
  enddo

  ! Bunit 

  b_unit_p(:)  = EXPRO_bunit(:)

  call EXPRO_alloc(path,0)

end subroutine cgyro_experimental_profiles
