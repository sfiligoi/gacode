!-------------------------------------------------------------------------
! qfm_interface.f90
!
! PURPOSE:
!  Provides interface description for QFM (Weiland model).
!
! CALLING SEQUENCE:
!  set qfm_*_in variables
!  call qfm_run(...)
!  get qfm_*_out variables
!-------------------------------------------------------------------------

module qfm_interface

  implicit none

  ! Input parameters

  real, dimension(27) :: qfm_in

  ! Output parameters

  real :: qfm_elec_pflux_out = 0.0
  real :: qfm_elec_eflux_out = 0.0
  real :: qfm_ion1_pflux_out = 0.0
  real :: qfm_ion1_eflux_out = 0.0
  real :: qfm_ion2_pflux_out = 0.0
  real :: qfm_ion2_eflux_out = 0.0
  real :: qfm_ion3_pflux_out = 0.0
  real :: qfm_ion3_eflux_out = 0.0

end module qfm_interface
