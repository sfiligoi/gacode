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

subroutine tgyro_expro_map(r,z,n_r,r_exp,p_exp,n_exp,mode)

  implicit none

  integer, intent(in) :: n_r,n_exp
  real, intent(in), dimension(n_r) :: r,z
  real, intent(in), dimension(n_exp) :: r_exp
  real, intent(inout), dimension(n_exp) :: p_exp
  real, dimension(n_exp) :: z_exp
  real :: dr
  integer :: i
  integer :: i_exp
  character(len=3) :: mode

  ! Compute z's on exp grid
  z_exp(:) = 0.0
  i = 1
  do i_exp=2,n_exp
     if (r_exp(i_exp) <= r(n_r)) then
        if (r(i+1) < r_exp(i_exp)) i = i+1
        dr = r(i+1)-r(i)
        z_exp(i_exp) = z(i+1)*(r_exp(i_exp)-r(i))/dr + z(i)*(r(i+1)-r_exp(i_exp))/dr
     endif
  enddo

  ! Start the integration at the first i_exp past r(n_r)
  ! ISSUE: This check is sensitive to mesh alignment (need =)
  do i_exp=n_exp,2,-1
     if (r_exp(i_exp-1) <= r(n_r)) exit
  enddo

  if (mode == 'log') then
     call logint(p_exp(1:i_exp),z_exp(1:i_exp),r_exp(1:i_exp),i_exp,i_exp)
  else
     call linint(p_exp(1:i_exp),z_exp(1:i_exp),r_exp(1:i_exp),i_exp,i_exp)
  endif

end subroutine tgyro_expro_map

