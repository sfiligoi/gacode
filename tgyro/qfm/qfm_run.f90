!---------------------------------------------------------
! qfm_run.f90
!---------------------------------------------------------

subroutine qfm_run

  use qfm_interface

  implicit none

  real, dimension(8) :: qfm_out

  call qfm_sub(qfm_in,qfm_out)

  qfm_elec_pflux_out = qfm_out(1)
  qfm_elec_eflux_out = qfm_out(2)
  qfm_ion1_pflux_out = qfm_out(3)
  qfm_ion1_eflux_out = qfm_out(4)
  qfm_ion2_pflux_out = qfm_out(5)
  qfm_ion2_eflux_out = qfm_out(6)
  qfm_ion3_pflux_out = qfm_out(7)
  qfm_ion3_eflux_out = qfm_out(8)

end subroutine qfm_run
