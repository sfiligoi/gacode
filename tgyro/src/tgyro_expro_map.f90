!-----------------------------------------------------------------
! tgyro_expro_map.f90
!
! PURPOSE:
!  Map evolved TGYRO profiles onto EXPRO grid *inside* r(n_r).
!
! NOTES:
!  Inputs:  r,z,n_r
!  Outputs: r_exp, p_exp, n_exp
!-----------------------------------------------------------------

subroutine tgyro_expro_map(r,z,n_r,ps,r_exp,p_exp,n_exp,mode)

  implicit none

  integer, intent(in) :: n_r,n_exp
  real, intent(in) :: ps
  real, intent(in), dimension(n_r) :: r,z
  real, intent(in), dimension(n_exp) :: r_exp
  real, intent(inout), dimension(n_exp) :: p_exp
  real, dimension(n_exp) :: z_exp
  real :: dr
  integer :: i
  integer :: i_exp
  character(len=3) :: mode

  do i_exp=1,n_exp
     if (r_exp(i_exp) <= r(n_r)) then
       call math_scaleint(z,r,n_r,ps,r_exp(i_exp),p_exp(i_exp),mode)
     endif
  enddo

end subroutine tgyro_expro_map

