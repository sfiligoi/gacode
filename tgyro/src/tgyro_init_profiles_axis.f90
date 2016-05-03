!----------------------------------------------------------------------
! tgyro_init_profiles_axis.f90
!
! PURPOSE:
!  Set profile values at r=0
!----------------------------------------------------------------------

subroutine tgyro_init_profiles_axis

  use tgyro_globals

  implicit none

  ! Zero derivatives of (ne,ni,Te,Ti) at origin.
  dlnnidr(:,1) = 0.0
  dlnnedr(1) = 0.0
  dlntidr(:,1) = 0.0
  dlntedr(1) = 0.0

  if (loc_er_feedback_flag == 1) then
     select case(tgyro_er_bc)
     case (1)
        ! Fixed derivative of w0 
        f_rot(1) = f_rot(1)
     case (2)
        ! Zero derivative of w0, similar to profiles
        f_rot(1) = 0.0
     case(3)
        ! Zero second derivative
        f_rot(1) = f_rot(2)
     end select
  endif

  pflux_e_neo(1) = 0.0
  pflux_e_tur(1) = 0.0
  pflux_i_neo(:,1) = 0.0
  pflux_i_tur(:,1) = 0.0

  eflux_e_neo(1) = 0.0
  eflux_e_tur(1) = 0.0
  eflux_i_neo(:,1) = 0.0
  eflux_i_tur(:,1) = 0.0

  mflux_e_neo(1) = 0.0
  mflux_e_tur(1) = 0.0
  mflux_i_neo(:,1) = 0.0
  mflux_i_tur(:,1) = 0.0

  pflux_i_tot(1) = 0.0
  pflux_e_tot(1) = 0.0

  eflux_i_tot(1) = 0.0
  eflux_e_tot(1) = 0.0

  mflux_tot(1) = 0.0
  pflux_he_tot(1) = 0.0

  eflux_i_target(1) = 0.0
  eflux_e_target(1) = 0.0
  pflux_e_target(1) = 0.0
  mflux_target(1)   = 0.0
  pflux_he_target(1) = 0.0

  ! Also need to zero initial exchanges to prevent use in tgyro_source 
  ! on iteration 0 before definition
  expwd_i_tur(:,:) = 0.0
  expwd_e_tur(:) = 0.0
  !----------------------------------------------------

end subroutine tgyro_init_profiles_axis
